from direct.showbase.ShowBase import ShowBase
from panda3d.core import *
from direct.gui.OnscreenText import OnscreenText
from direct.task import Task
import math
import numpy as np
from scipy.optimize import minimize
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C


class ScipyOptimizer:
    """
    Nelder-Mead (Downhill Simplex) Optimization for 3D laser alignment
    
    Uses scipy.optimize.minimize with method='Nelder-Mead' to find optimal mirror angles.
    Wraps the simulation's intensity evaluation function to provide scipy with a callable objective.
    
    Key Features:
    - No gradient computation needed - works with noisy objective functions
    - Maintains simplex of n+1 points in n-dimensional space (5 points for 4D)
    - Adapts to local topology through reflection, expansion, contraction operations
    - Typically converges in 50-200 evaluations for this problem
    - Captures full optimization history for animation/analysis
    
    Advantages over other methods:
    - Faster than grid search (100-200 vs 10,000+ samples)
    - More robust than Powell for non-separable functions
    - Simpler than Bayesian optimization (no GP fitting)
    - Works well for quasi-convex landscapes like fiber coupling
    """
    
    def __init__(self, maxiter=200, bounds_deg=15, n_restarts=3, shrink_factor=0.1):
        """
        Args:
            maxiter: Maximum number of iterations per restart (function evaluations)
            bounds_deg: Maximum rotation from initial position in degrees
            n_restarts: Number of restarts with shrinking simplexes
            shrink_factor: Factor to shrink simplex between restarts (e.g., 0.5 = half size)
        """
        self.maxiter = maxiter
        self.bounds_deg = bounds_deg
        self.n_restarts = n_restarts
        self.shrink_factor = shrink_factor
        
        self.is_running = False
        self.best_intensity = -1
        self.best_orientations = None
        self.best_angles = None
        
        # Optimization history
        self.history = []  # List of {'iteration': int, 'angles': array, 'intensity': float}
        self.current_iteration = 0
        
        # Current optimization state
        self.initial_quat1 = None
        self.initial_quat2 = None
        
        # Evaluation callback (set by simulation)
        self.evaluate_intensity = None
        
        # Scipy optimization result
        self.scipy_result = None
        self.optimization_complete = False
        
        # Reference to simulation app (for processing graphics during optimization)
        self.app = None
        
    def start(self, disk1, disk2, evaluate_callback, app=None):
        """Start optimization from current disk orientations
        
        Args:
            disk1: First mirror disk node
            disk2: Second mirror disk node
            evaluate_callback: Function that takes angles [ax1, ay1, ax2, ay2] and returns intensity
            app: Reference to Panda3D ShowBase app for processing graphics updates
        """
        self.is_running = True
        self.initial_quat1 = disk1.getQuat()
        self.initial_quat2 = disk2.getQuat()
        self.best_orientations = (self.initial_quat1, self.initial_quat2)
        self.best_intensity = -1
        self.evaluate_intensity = evaluate_callback
        self.app = app
        
        # Reset state
        self.history = []
        self.current_iteration = 0
        self.optimization_complete = False
        self.last_angles = None
        
        # Run scipy optimization (this will call evaluate_callback many times)
        self._run_scipy_optimization()
        
    def stop(self):
        """Stop optimization"""
        self.is_running = False
        
    def _run_scipy_optimization(self):
        """Run scipy.optimize.minimize with Nelder-Mead method with multiple restarts"""
        print(f"NelderMeadOptimizer: Starting scipy.optimize.minimize with {self.n_restarts} restarts...")
        
        # Bounds
        max_rad = math.radians(self.bounds_deg)
        
        # Objective function: minimize negative intensity (to maximize intensity)
        def objective(angles):
            # Evaluate intensity through simulation
            intensity = self.evaluate_intensity(angles)
            
            # Process graphics updates so we can see the optimization happening
            # Render one frame without processing tasks (avoids recursive task warnings)
            if self.app is not None:
                self.app.graphicsEngine.renderFrame()
            
            actual_angles = getattr(self, 'last_angles', np.array(angles, dtype=float))

            # Store in history
            self.history.append({
                'iteration': self.current_iteration,
                'angles': actual_angles.copy(),
                'intensity': intensity,
                'alpha1': actual_angles[0],
                'theta1': actual_angles[1],
                'alpha2': actual_angles[2],
                'theta2': actual_angles[3]
            })
            # Debug output every 25 iterations
            if self.current_iteration % 25 == 0:
                print(f"  Iter {self.current_iteration}: intensity = {intensity:.4f}")
            
            self.current_iteration += 1
            
            # Update best
            if intensity > self.best_intensity:
                self.best_intensity = intensity
                self.best_orientations = self._angles_to_quaternions(actual_angles)
                self.best_angles = actual_angles.copy()
            
            # Return negative for minimization
            return 1-intensity
        
        # Multi-restart optimization with exponentially shrinking simplexes
        for restart in range(self.n_restarts):
            # Calculate simplex scale for this restart
            simplex_scale = max_rad * (self.shrink_factor ** restart)
            
            # Initial guess: start from best found so far (or zero for first restart)
            if restart == 0:
                initial_guess = np.zeros(4)
            else:
                # Start from best angles found so far
                initial_guess = self.best_angles.copy() if self.best_angles is not None else np.zeros(4)
            
            self.last_angles = initial_guess.copy()
            
            # Configure initial simplex for this restart
            n = len(initial_guess)
            p1 = -1.0/(n*math.sqrt(2)) * (math.sqrt(n+1) - 1)
            p2 = 1.0/(n*math.sqrt(2)) * (math.sqrt(n+1) + n - 1)
            
            # Create initial simplex with current scale
            initial_simplex = np.zeros((n + 1, n))
            initial_simplex[0] = initial_guess + np.full(n, p1) * simplex_scale
            for i in range(n):
                modifier = 2 if i < 2 else 1
                point = initial_guess + np.full(n, p1) * simplex_scale
                point[i] = initial_guess[i] + p2 * simplex_scale * modifier
                initial_simplex[i + 1] = point
            
            print(f"\n  Restart {restart + 1}/{self.n_restarts}: simplex_scale = {math.degrees(simplex_scale):.4f}°")
            print(f"    Starting from: {np.rad2deg(initial_guess)} deg")
            
            # Run optimization for this restart
            result = minimize(
                objective,
                initial_guess,
                method='Nelder-Mead',
                options={
                    'maxfev': self.maxiter,
                    # 'xatol': math.radians(2*np.pi/2048 / 500),
                    'xatol': 1e-100,
                    'fatol': 1e-100,
                    'initial_simplex': initial_simplex,
                    'adaptive': False
                }
            )
            
            print(f"    Restart {restart + 1} complete: best_intensity = {self.best_intensity:.4f}")
            
            # Store the final result from last restart
            if restart == self.n_restarts - 1:
                self.scipy_result = result
        
        print(f"\nNelderMeadOptimizer: All restarts complete!")
        print(f"  Total iterations: {len(self.history)}")
        print(f"  Best intensity: {self.best_intensity:.4f}")
        if self.best_angles is not None:
            print(f"  Best angles: {np.rad2deg(self.best_angles)} deg")
        
        self.optimization_complete = True
        
    def _angles_to_quaternions(self, angles):
        """Convert angle parameters to quaternions"""
        ax1, ay1, ax2, ay2 = angles
        
        # Create rotation quaternions for disk 1
        rot1_x = LQuaternionf()
        rot1_x.setFromAxisAngle(math.degrees(ax1), LVector3f(1, 0, 0))
        rot1_y = LQuaternionf()
        rot1_y.setFromAxisAngle(math.degrees(ay1), LVector3f(0, 1, 0))
        
        # Create rotation quaternions for disk 2
        rot2_x = LQuaternionf()
        rot2_x.setFromAxisAngle(math.degrees(ax2), LVector3f(1, 0, 0))
        rot2_y = LQuaternionf()
        rot2_y.setFromAxisAngle(math.degrees(ay2), LVector3f(0, 1, 0))
        
        # Combine with initial orientations
        quat1 = rot1_y * rot1_x * self.initial_quat1
        quat2 = rot2_y * rot2_x * self.initial_quat2
        
        return quat1, quat2
    
    def get_best_state(self):
        """Get best orientation found so far"""
        if len(self.history) == 0:
            return None
        
        # Find best in history
        best_entry = max(self.history, key=lambda h: h['intensity'])
        best_angles = best_entry['angles']
        
        return self._angles_to_quaternions(best_angles)
        
    def is_complete(self):
        """Check if optimization is complete"""
        return self.optimization_complete or not self.is_running


class GaussianProcessOptimizer:
    """
    Gaussian Process (Bayesian) Optimization for 3D laser alignment
    
    Uses a Gaussian Process to model the intensity function and intelligently
    select the next points to sample based on Expected Improvement (EI).
    
    Strategy:
    1. Random exploration: Sample random points to build initial GP model
    2. Exploitation: Use GP to find points with high expected improvement
    3. Balance exploration/exploitation through acquisition function
    
    Advantages:
    - Very sample-efficient (often finds optimum in 30-100 samples)
    - Models uncertainty, not just the mean prediction
    - Good for expensive black-box functions
    - Adapts to local landscape automatically
    
    How it works:
    - GP models intensity as f(angles) ~ N(μ(angles), σ²(angles))
    - High μ = high predicted intensity
    - High σ = high uncertainty (unexplored region)
    - EI balances both: try high μ OR high σ regions
    """
    
    def __init__(self, n_initial_samples=30, n_total_samples=100, bounds_deg=15):
        """
        Args:
            n_initial_samples: Random samples before using GP
            n_total_samples: Total samples to take
            bounds_deg: Maximum rotation from initial position in degrees
        """
        self.n_initial_samples = n_initial_samples
        self.n_total_samples = n_total_samples
        self.bounds_deg = bounds_deg
        
        self.is_running = False
        self.best_intensity = -1
        self.best_orientations = None
        
        # GP model
        self.gp = None
        self.X_samples = []  # List of sampled angles
        self.y_samples = []  # List of observed intensities
        
        # Optimization state
        self.current_sample = 0
        self.initial_quat1 = None
        self.initial_quat2 = None
        self.evaluate_intensity = None
        self.app = None
        
        # For async operation
        self.current_angles = None
        self.optimization_complete = False
        
    def start(self, disk1, disk2, evaluate_callback, app=None):
        """Start optimization from current disk orientations"""
        self.is_running = True
        self.initial_quat1 = disk1.getQuat()
        self.initial_quat2 = disk2.getQuat()
        self.best_orientations = (self.initial_quat1, self.initial_quat2)
        self.best_intensity = -1
        self.evaluate_intensity = evaluate_callback
        self.app = app
        
        # Reset state
        self.X_samples = []
        self.y_samples = []
        self.current_sample = 0
        self.optimization_complete = False
        
        # Initialize GP with RBF kernel optimized for smooth, unimodal functions
        # Longer length scale = smoother function assumption = more convex-like
        kernel = C(1.0, (0.1, 10.0)) * RBF(length_scale=5.0, length_scale_bounds=(1e-5, 1e5))
        self.gp = GaussianProcessRegressor(
            kernel=kernel,
            n_restarts_optimizer=10,  # More restarts for kernel hyperparameter optimization
            alpha=1e-6,
            normalize_y=True
        )
        
        print(f"GaussianProcessOptimizer: Starting with {self.n_initial_samples} random samples")
        
    def stop(self):
        """Stop optimization"""
        self.is_running = False
        
    def get_next_state(self):
        """Get next state to try (called each frame)"""
        if not self.is_running or self.optimization_complete:
            return None
            
        if self.current_sample >= self.n_total_samples:
            self.optimization_complete = True
            print(f"\nGaussianProcessOptimizer: Complete!")
            print(f"  Total samples: {len(self.X_samples)}")
            print(f"  Best intensity: {self.best_intensity:.4f}")
            if len(self.X_samples) > 0:
                best_idx = np.argmax(self.y_samples)
                best_angles = self.X_samples[best_idx]
                print(f"  Best angles: {np.rad2deg(best_angles)} deg")
                self.best_orientations = self._angles_to_quaternions(best_angles)
            return None
        
        # Decide next point to sample
        # if self.current_sample < self.n_initial_samples:
        #     # Random exploration phase
        #     angles = self._sample_random_point()
        # else:
        #     # GP-guided exploitation phase
        #     angles = self._sample_next_point_ei()
        angles = self._find_best_predicted_point()
        
        # Evaluate intensity
        intensity = self.evaluate_intensity(angles)
        
        # Store sample (store all samples)
        if intensity > 0:  # ignore zero intensities because they represent laser missing sensor
            self.X_samples.append(angles)
            self.y_samples.append(intensity)
            self.current_sample += 1
        
        # Update best
        if intensity > self.best_intensity:
            self.best_intensity = intensity
            self.best_orientations = self._angles_to_quaternions(angles)
        
        # Update GP model (after we have enough samples)
        if self.current_sample >= self.n_initial_samples and len(self.X_samples) >= self.n_initial_samples:
            X = np.array(self.X_samples)
            y = np.array(self.y_samples)
            try:
                self.gp.fit(X, y)
                
                # Debug: Check Gaussian fit quality
                if self.current_sample % 10 == 0:
                    # Predict on training data to check fit quality
                    y_pred, y_std = self.gp.predict(X, return_std=True)
                    
                    # Calculate R² score (how well GP fits the data)
                    residuals = y - y_pred
                    ss_res = np.sum(residuals**2)
                    ss_tot = np.sum((y - np.mean(y))**2)
                    r2_score = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
                    
                    # Log marginal likelihood (higher is better fit)
                    log_likelihood = self.gp.log_marginal_likelihood_value_
                    
                    # Mean prediction uncertainty
                    mean_std = np.mean(y_std)
                    
                    # Check residual normality: should be roughly N(0, σ²)
                    mean_residual = np.mean(residuals)
                    std_residual = np.std(residuals)
                    
                    print(f"  [GP Fit Quality] R²={r2_score:.4f} | LogML={log_likelihood:.2f} | "
                          f"MeanStd={mean_std:.4f} | ResidMean={mean_residual:.4f}±{std_residual:.4f}")
                    
            except Exception as e:
                print(f"  Warning: GP fit failed: {e}")
        
        # Render frame
        if self.app is not None:
            self.app.graphicsEngine.renderFrame()
        
        if self.current_sample % 10 == 0:
            phase = "Explore" if self.current_sample < self.n_initial_samples else "Exploit"
            print(f"  GP {phase} sample {self.current_sample}: intensity = {intensity:.4f}, best = {self.best_intensity:.4f}")
        
        
        # Return quaternions for this sample
        return self._angles_to_quaternions(angles)
    
    def _sample_random_point(self):
        """Sample a random point within bounds"""
        max_rad = math.radians(self.bounds_deg)
        angles = np.random.uniform(-max_rad, max_rad, size=4)
        return angles
    
    def _find_best_predicted_point(self):
        """Find the angle that maximizes the GP's predicted mean intensity
        
        Uses scipy.optimize.minimize to find the maximum of the GP's mean prediction.
        Since the GP fit is excellent and assumes a smooth/convex-like function,
        a single optimization from the best observed point should find the global maximum.
        
        Returns:
            Array of angles [ax1, ay1, ax2, ay2] that maximize predicted intensity
        """
        if self.gp is None or len(self.X_samples) < self.n_initial_samples:
            # Fall back to random sampling if GP not ready
            return self._sample_random_point()
        
        max_rad = math.radians(self.bounds_deg)
        
        # Objective: minimize negative predicted mean (to maximize mean)
        def objective(angles):
            angles_2d = angles.reshape(1, -1)
            mu = self.gp.predict(angles_2d, return_std=False)
            return -mu[0]  # Negative for minimization
        
        # Start from best observed sample
        best_idx = np.argmax(self.y_samples)
        initial_guess = np.array(self.X_samples[best_idx])
        
        # Single optimization run (no multi-restart needed with smooth GP)
        result = minimize(
            objective,
            initial_guess,
            method='L-BFGS-B',
            bounds=[(-max_rad, max_rad)] * 4,
            options={
                'maxiter': 100,
                'ftol': 1e-9
            }
        )
        
        return result.x
    def _angles_to_quaternions(self, angles):
        """Convert angle parameters to quaternions"""
        ax1, ay1, ax2, ay2 = angles
        
        # Create rotation quaternions for disk 1
        rot1_x = LQuaternionf()
        rot1_x.setFromAxisAngle(math.degrees(ax1), LVector3f(1, 0, 0))
        rot1_y = LQuaternionf()
        rot1_y.setFromAxisAngle(math.degrees(ay1), LVector3f(0, 1, 0))
        
        # Create rotation quaternions for disk 2
        rot2_x = LQuaternionf()
        rot2_x.setFromAxisAngle(math.degrees(ax2), LVector3f(1, 0, 0))
        rot2_y = LQuaternionf()
        rot2_y.setFromAxisAngle(math.degrees(ay2), LVector3f(0, 1, 0))
        
        # Combine with initial orientations
        quat1 = rot1_y * rot1_x * self.initial_quat1
        quat2 = rot2_y * rot2_x * self.initial_quat2
        
        return quat1, quat2
    
    def get_best_state(self):
        """Get best orientation found so far"""
        return self.best_orientations
        
    def is_complete(self):
        """Check if optimization is complete"""
        return self.optimization_complete or not self.is_running


class LaserSimulation(ShowBase):
    def __init__(self):
        ShowBase.__init__(self)
        
        # Window setup
        self.win.setClearColor((0.1, 0.1, 0.15, 1))
        
        # Setup camera
        self.disableMouse()
        
        # Camera control variables
        self.camera_distance = 30
        self.camera_heading = 0  # Horizontal rotation (around Z axis)
        self.camera_pitch = 30   # Vertical rotation
        self.camera_target = LVector3f(0, 0, 0)
        
        # Mouse control state
        self.is_orbiting = False
        self.is_panning = False
        self.last_mouse_x = 0
        self.last_mouse_y = 0
        
        # Update camera position
        self.update_camera_position()
        
        # Setup lighting
        self.setup_lights()
        
        # Initialize scene
        wire_pos = LVector3f(-50, 0, 0)
        laser_pos = LVector3f(0, -10, 0)
        self.laser_pos = LVector3f(0, -10, 0)
        self.laser_dir = LVector3f(1, 0, 0)
        self.laser_dir.normalize()
        
        # Create disk mirrors
        disk1_center = LVector3f(10, -10, 0)
        disk2_center = LVector3f(12, 1, 0)
        
        # Mirror 1 normal: pointing toward average of laser source and mirror 2
        disk1_normal = (laser_pos + disk2_center)/2 - disk1_center
        disk1_normal.normalize()
        
        # Mirror 2 normal: bisector of incident ray from mirror 1 and direction to wire
        # First get the reflected direction from mirror 1
        to_disk2 = disk2_center - disk1_center
        to_disk2.normalize()
        from_disk1 = -to_disk2  # Approximate: ray coming from disk1
        to_wire = wire_pos - disk2_center
        to_wire.normalize()
        disk2_normal = from_disk1 + to_wire
        disk2_normal.normalize()

        self.disk1 = self.create_disk_mirror(
            center=disk1_center,
            normal=disk1_normal,
            radius=4.0,
            color=(0.8, 0.8, 1, 1)
        )
        self.disk2 = self.create_disk_mirror(
            center=disk2_center,
            normal=disk2_normal,
            radius=4.0,
            color=(0.8, 0.8, 1, 1)
        )
        
        self.disks = [self.disk1, self.disk2]
        
        # Store initial orientations for clamping
        self.disk1_initial_quat = self.disk1.getQuat()
        self.disk2_initial_quat = self.disk2.getQuat()
        self.max_rotation_deg = 10  # Maximum rotation from initial position
        
        # Create wire detector
        self.wire = self.create_wire_detector(
            position=wire_pos,
            direction=LVector3f(1, 0, 0),
            radius=0.5,
            length=3.0
        )
        
        # Create ground plane
        self.create_ground()
        
        # Ray path lines
        self.ray_lines = None
        self.ray_line_node = None
        
        # UI Text
        self.intensity_text = OnscreenText(
            text='Wire Coupling Intensity: 0.0000',
            pos=(-1.3, 0.9),
            scale=0.05,
            fg=(1, 1, 0, 1),
            align=TextNode.ALeft
        )  
        self.optimizer_text = OnscreenText(
            text='Optimizer: OFF',
            pos=(-1.3, 0.8),
            scale=0.04,
            fg=(0, 1, 1, 1),
            align=TextNode.ALeft
        )  
        self.controls_text = OnscreenText(
            text='Disk 1: Arrow Keys | Disk 2: WASD | Camera: Left-Click Orbit, Scroll Zoom, Middle-Click Pan | O: Optimize',
            pos=(-1.3, -0.9),
            scale=0.04,
            fg=(1, 1, 1, 1),
            align=TextNode.ALeft
        )
        
        # Keyboard state
        self.keys = {}
        self.setup_controls()
    
        # Choose optimizer:
        # Full 4D optimization (may take 200-500 evaluations):
        self.optimizer = ScipyOptimizer(maxiter=500, bounds_deg=3, n_restarts=1, shrink_factor=0.3)
        
        
        # Gaussian Process optimization (typically 30-100 evaluations):
        # Uses Bayesian optimization with Expected Improvement
        # self.optimizer = GaussianProcessOptimizer(n_initial_samples=15, n_total_samples=150, bounds_deg=1/10)
        
        # Add update task
        self.taskMgr.add(self.update_task, "update")
        
    def update_camera_position(self):
        """Update camera position based on orbit parameters"""
        # Convert spherical coordinates to cartesian
        heading_rad = math.radians(self.camera_heading)
        pitch_rad = math.radians(self.camera_pitch)
        
        # Calculate camera position
        x = self.camera_distance * math.cos(pitch_rad) * math.sin(heading_rad)
        y = self.camera_distance * math.cos(pitch_rad) * math.cos(heading_rad)
        z = self.camera_distance * math.sin(pitch_rad)
        
        self.camera.setPos(
            self.camera_target.x + x,
            self.camera_target.y + y,
            self.camera_target.z + z
        )
        self.camera.lookAt(self.camera_target)
    
    def setup_lights(self):
        """Setup scene lighting"""
        # Ambient light
        alight = AmbientLight('alight')
        alight.setColor((0.3, 0.3, 0.3, 1))
        alnp = self.render.attachNewNode(alight)
        self.render.setLight(alnp)
        
        # Directional light 1
        dlight1 = DirectionalLight('dlight1')
        dlight1.setColor((0.8, 0.8, 0.8, 1))
        dlnp1 = self.render.attachNewNode(dlight1)
        dlnp1.setHpr(45, -45, 0)
        self.render.setLight(dlnp1)
        
        # Directional light 2
        dlight2 = DirectionalLight('dlight2')
        dlight2.setColor((0.5, 0.5, 0.5, 1))
        dlnp2 = self.render.attachNewNode(dlight2)
        dlnp2.setHpr(-45, -30, 0)
        self.render.setLight(dlnp2)
        
    def create_disk_mirror(self, center, normal, radius, color):
        """Create a circular mirror disk"""
        # Create circle mesh
        format = GeomVertexFormat.getV3n3c4()
        vdata = GeomVertexData('disk', format, Geom.UHStatic)
        vdata.setNumRows(33)
        
        vertex = GeomVertexWriter(vdata, 'vertex')
        normal_writer = GeomVertexWriter(vdata, 'normal')
        color_writer = GeomVertexWriter(vdata, 'color')
        
        # Center vertex
        vertex.addData3(0, 0, 0)
        normal_writer.addData3(0, 0, 1)
        color_writer.addData4(*color)
        
        # Circle vertices
        segments = 32
        for i in range(segments + 1):
            angle = 2 * math.pi * i / segments
            x = radius * math.cos(angle)
            y = radius * math.sin(angle)
            vertex.addData3(x, y, 0)
            normal_writer.addData3(0, 0, 1)
            color_writer.addData4(*color)
        
        # Create triangles
        tris = GeomTriangles(Geom.UHStatic)
        for i in range(segments):
            tris.addVertices(0, i + 1, i + 2)
        
        geom = Geom(vdata)
        geom.addPrimitive(tris)
        
        node = GeomNode('disk')
        node.addGeom(geom)
        
        disk_np = self.render.attachNewNode(node)
        disk_np.setPos(center)
        disk_np.setTwoSided(True)
        
        # Orient disk to face normal direction
        # The disk is created in XY plane with normal along +Z
        # We need to rotate it so its local +Z points along 'normal'
        normal.normalize()
        
        # Use lookAt to orient the disk - the disk's +Z should point along normal
        # lookAt makes the node's +Y axis point at the target by default
        # So we need to rotate 90 degrees to make +Z point at the target
        disk_np.lookAt(center + normal)
        disk_np.setP(disk_np.getP() - 90)  # Rotate to make +Z (not +Y) point at target
        
        # Store properties
        disk_np.setPythonTag('normal', normal)
        disk_np.setPythonTag('radius', radius)
        disk_np.setPythonTag('center', center)
        
        # Add normal indicator
        # self.create_normal_indicator(disk_np, center, normal, radius)
        
        return disk_np
    
    def create_normal_indicator(self, parent, center, normal, radius):
        """Create a visual indicator for the disk normal"""
        return
        # Create a small cylinder pointing along normal
        from panda3d.core import CardMaker
        cm = CardMaker('normal_indicator')
        cm.setFrame(-0.1, 0.1, 0, radius * 0.6)
        indicator = parent.attachNewNode(cm.generate())
        indicator.setColor(0, 1, 1, 1)  # Cyan
        indicator.setP(-90)  # Point forward
        
    def create_wire_detector(self, position, direction, radius, length):
        """Create wire detector (fiber optic)"""
        # Create a box for the wire
        format = GeomVertexFormat.getV3n3c4()
        vdata = GeomVertexData('wire', format, Geom.UHStatic)
        
        vertex = GeomVertexWriter(vdata, 'vertex')
        normal_writer = GeomVertexWriter(vdata, 'normal')
        color_writer = GeomVertexWriter(vdata, 'color')
        
        # Simple box vertices
        hw, hh, hl = length / 2, radius, radius
        vertices = [
            (-hl, -hw, -hh), (hl, -hw, -hh), (hl, hw, -hh), (-hl, hw, -hh),  # Back
            (-hl, -hw, hh), (hl, -hw, hh), (hl, hw, hh), (-hl, hw, hh),  # Front
        ]
        
        for v in vertices:
            vertex.addData3(*v)
            normal_writer.addData3(0, 0, 1)
            color_writer.addData4(0.2, 0.8, 0.2, 1)
        
        # Create box faces
        tris = GeomTriangles(Geom.UHStatic)
        faces = [
            (0, 1, 2), (0, 2, 3),  # Back
            (4, 5, 6), (4, 6, 7),  # Front
            (0, 1, 5), (0, 5, 4),  # Bottom
            (2, 3, 7), (2, 7, 6),  # Top
            (0, 3, 7), (0, 7, 4),  # Left
            (1, 2, 6), (1, 6, 5),  # Right
        ]
        for face in faces:
            tris.addVertices(*face)
        
        geom = Geom(vdata)
        geom.addPrimitive(tris)
        
        node = GeomNode('wire')
        node.addGeom(geom)
        
        wire_np = self.render.attachNewNode(node)
        wire_np.setPos(position)
        wire_np.setTwoSided(True)
        
        direction.normalize()
        wire_np.lookAt(position + direction)
        wire_np.setR(90)
        
        # Store properties
        wire_np.setPythonTag('position', position)
        wire_np.setPythonTag('direction', direction)
        wire_np.setPythonTag('radius', radius)
        wire_np.setPythonTag('length', length)
        
        # Create entrance indicator (yellow circle)
        self.create_entrance_indicator(wire_np, position, radius)
        
        return wire_np
    
    def create_entrance_indicator(self, parent, position, radius):
        """Create entrance circle indicator"""
        format = GeomVertexFormat.getV3n3c4()
        vdata = GeomVertexData('entrance', format, Geom.UHStatic)
        
        vertex = GeomVertexWriter(vdata, 'vertex')
        normal_writer = GeomVertexWriter(vdata, 'normal')
        color_writer = GeomVertexWriter(vdata, 'color')
        
        # Center
        vertex.addData3(0, 0, 0)
        normal_writer.addData3(0, 0, 1)
        color_writer.addData4(1, 1, 0, 1)
        
        # Circle
        segments = 16
        for i in range(segments + 1):
            angle = 2 * math.pi * i / segments
            x = radius * math.cos(angle)
            y = radius * math.sin(angle)
            vertex.addData3(0, x, y)
            normal_writer.addData3(1, 0, 0)
            color_writer.addData4(1, 1, 0, 1)
        
        tris = GeomTriangles(Geom.UHStatic)
        for i in range(segments):
            tris.addVertices(0, i + 1, i + 2)
        
        geom = Geom(vdata)
        geom.addPrimitive(tris)
        
        node = GeomNode('entrance')
        node.addGeom(geom)
        
        entrance_np = parent.attachNewNode(node)
        entrance_np.setPos(-radius * 1.5, 0, 0)
        
    def create_ground(self):
        """Create ground plane"""
        from panda3d.core import CardMaker
        cm = CardMaker('ground')
        cm.setFrame(-2500, 2500, -2500, 2500)
        ground = self.render.attachNewNode(cm.generate())
        ground.setP(-90)
        ground.setPos(0, 0, -5)
        ground.setColor(0.2, 0.2, 0.2, 1)
        
    def setup_controls(self):
        """Setup keyboard controls"""
        self.accept('arrow_left', self.set_key, ['left', True])
        self.accept('arrow_left-up', self.set_key, ['left', False])
        self.accept('arrow_right', self.set_key, ['right', True])
        self.accept('arrow_right-up', self.set_key, ['right', False])
        self.accept('arrow_up', self.set_key, ['up', True])
        self.accept('arrow_up-up', self.set_key, ['up', False])
        self.accept('arrow_down', self.set_key, ['down', True])
        self.accept('arrow_down-up', self.set_key, ['down', False])
        
        self.accept('w', self.set_key, ['w', True])
        self.accept('w-up', self.set_key, ['w', False])
        self.accept('s', self.set_key, ['s', True])
        self.accept('s-up', self.set_key, ['s', False])
        self.accept('a', self.set_key, ['a', True])
        self.accept('a-up', self.set_key, ['a', False])
        self.accept('d', self.set_key, ['d', True])
        self.accept('d-up', self.set_key, ['d', False])
        
        self.accept('o', self.toggle_optimizer)
        
        # Mouse controls
        self.accept('mouse1', self.start_orbit)
        self.accept('mouse1-up', self.stop_orbit)
        self.accept('mouse2', self.start_pan)
        self.accept('mouse2-up', self.stop_pan)
        self.accept('wheel_up', self.zoom_in)
        self.accept('wheel_down', self.zoom_out)
        
    def set_key(self, key, value):
        """Set keyboard state"""
        self.keys[key] = value
        
    def start_orbit(self):
        """Start camera orbit on left mouse button"""
        self.is_orbiting = True
        if self.mouseWatcherNode.hasMouse():
            self.last_mouse_x = self.mouseWatcherNode.getMouseX()
            self.last_mouse_y = self.mouseWatcherNode.getMouseY()
            
    def stop_orbit(self):
        """Stop camera orbit"""
        self.is_orbiting = False
        
    def start_pan(self):
        """Start camera panning on middle mouse button"""
        self.is_panning = True
        if self.mouseWatcherNode.hasMouse():
            self.last_mouse_x = self.mouseWatcherNode.getMouseX()
            self.last_mouse_y = self.mouseWatcherNode.getMouseY()
            
    def stop_pan(self):
        """Stop camera panning"""
        self.is_panning = False
        
    def zoom_in(self):
        """Zoom camera in"""
        self.camera_distance = max(5, self.camera_distance - 2)
        self.update_camera_position()
        
    def zoom_out(self):
        """Zoom camera out"""
        self.camera_distance = min(100, self.camera_distance + 2)
        self.update_camera_position()
        
    def toggle_optimizer(self):
        """Toggle the optimizer on/off"""
        if self.optimizer.is_running:
            self.optimizer.stop()
            self.optimizer_text.setText('Optimizer: OFF')
        else:
            # Start optimizer - pass evaluation callback for optimizers that need it
            if isinstance(self.optimizer, ScipyOptimizer):
                self.optimizer.start(self.disk1, self.disk2, self.evaluate_intensity_for_angles, app=self)
                opt_name = "Nelder-Mead"
            elif isinstance(self.optimizer, GaussianProcessOptimizer):
                self.optimizer.start(self.disk1, self.disk2, self.evaluate_intensity_for_angles, app=self)
                opt_name = "Gaussian Process"
            elif isinstance(self.optimizer, BayesianOptimizer):
                self.optimizer.start(self.disk1, self.disk2)
                opt_name = "Bayesian"
            elif isinstance(self.optimizer, CoordinateGridSearch):
                self.optimizer.start(self.disk1, self.disk2)
                opt_name = "Coordinate Grid Search"
            elif isinstance(self.optimizer, PowellOptimizer):
                self.optimizer.start(self.disk1, self.disk2)
                opt_name = "Powell"
            else:
                self.optimizer.start(self.disk1, self.disk2)
                opt_name = "Grid Search"
            self.optimizer_text.setText(f'{opt_name} Optimizer: STARTING...')
            
    def clamp_disk_rotation(self, disk, initial_quat):
        """Clamp disk rotation to max_rotation_deg from initial orientation"""
        return
        current_quat = disk.getQuat()
        
        # Calculate rotation difference
        # q_diff = q_initial^-1 * q_current
        initial_quat_inv = LQuaternionf()
        initial_quat_inv.invertFrom(initial_quat)
        diff_quat = initial_quat_inv * current_quat
        
        # Get angle of rotation
        angle = diff_quat.getAngle()
        
        # If rotation exceeds limit, clamp it
        if angle > self.max_rotation_deg:
            # Scale down the rotation
            axis = diff_quat.getAxis()
            axis.normalize()  # Ensure axis is normalized
            
            # Check if axis is valid
            if axis.lengthSquared() < 0.001:
                # If axis is too small, just use the initial orientation
                disk.setQuat(initial_quat)
            else:
                clamped_diff = LQuaternionf()
                clamped_diff.setFromAxisAngle(self.max_rotation_deg, axis)
                
                # Apply clamped rotation
                clamped_quat = initial_quat * clamped_diff
                disk.setQuat(clamped_quat)
            
            # Update normal vector
            forward = LVector3f(0, 0, 1)
            new_normal = disk.getQuat().xform(forward)
            new_normal.normalize()
            disk.setPythonTag('normal', new_normal)
            
            return True  # Rotation was clamped
        return False  # No clamping needed
        
    def rotate_disk(self, disk, axis, angle_deg):
        """Rotate a disk around an axis"""
        # Get current orientation
        current_quat = disk.getQuat()
        
        # Create rotation quaternion
        if axis == 'x':
            rot_quat = LQuaternionf()
            rot_quat.setFromAxisAngle(angle_deg, LVector3f(1, 0, 0))
        elif axis == 'y':
            rot_quat = LQuaternionf()
            rot_quat.setFromAxisAngle(angle_deg, LVector3f(0, 1, 0))
        else:
            rot_quat = LQuaternionf()
            rot_quat.setFromAxisAngle(angle_deg, LVector3f(0, 0, 1))
        
        # Apply rotation
        new_quat = current_quat * rot_quat
        disk.setQuat(new_quat)
        
        # Clamp rotation
        initial_quat = self.disk1_initial_quat if disk == self.disk1 else self.disk2_initial_quat
        self.clamp_disk_rotation(disk, initial_quat)
        
        # Update normal vector
        forward = LVector3f(0, 0, 1)
        new_normal = disk.getQuat().xform(forward)
        new_normal.normalize()
        disk.setPythonTag('normal', new_normal)
        
    def calculate_disk_reflection(self, ray_pos, ray_dir, disk):
        """Calculate ray-disk intersection and reflection"""
        n = disk.getPythonTag('normal')
        c = disk.getPythonTag('center')
        radius = disk.getPythonTag('radius')
        
        # Update center to current position
        c = disk.getPos()
        
        to_center = c - ray_pos
        denom = ray_dir.dot(n)
        
        if abs(denom) < 0.001:
            return None
        
        t = to_center.dot(n) / denom
        if t < 0.01:
            return None
        
        intersection = ray_pos + ray_dir * t
        to_intersection = intersection - c
        dist_sq = to_intersection.lengthSquared()
        
        if dist_sq > radius ** 2:
            return None
        
        # Calculate reflection
        dot_val = ray_dir.dot(n)
        reflection = ray_dir - n * (2 * dot_val)
        reflection.normalize()
        
        return intersection, reflection, math.sqrt(dist_sq)
    
    def trace_ray_3d(self):
        """Trace the laser ray through the scene"""
        current_pos = LVector3f(self.laser_pos)
        current_dir = LVector3f(self.laser_dir)
        path = [LVector3f(current_pos)]
        max_bounces = 10
        
        for _ in range(max_bounces):
            closest_dist = float('inf')
            closest_result = None
            
            for disk in self.disks:
                result = self.calculate_disk_reflection(current_pos, current_dir, disk)
                if result:
                    intersection, reflection, _ = result
                    dist = (intersection - current_pos).length()
                    if dist < closest_dist:
                        closest_dist = dist
                        closest_result = (intersection, reflection)
            
            if closest_result:
                intersection, reflection = closest_result
                path.append(LVector3f(intersection))
                current_pos = intersection
                current_dir = reflection
            else:
                break
        
        # Extend ray beyond last bounce
        end_pos = current_pos + current_dir * 2000
        path.append(end_pos)
        
        return path, current_pos, current_dir
    
    def draw_ray_path(self, path):
        """Draw the ray path as lines"""
        if self.ray_line_node:
            self.ray_line_node.removeNode()
        
        if len(path) < 2:
            return
        
        # Create line segments
        format = GeomVertexFormat.getV3c4()
        vdata = GeomVertexData('ray', format, Geom.UHDynamic)
        
        vertex = GeomVertexWriter(vdata, 'vertex')
        color_writer = GeomVertexWriter(vdata, 'color')
        
        for point in path:
            vertex.addData3(point)
            color_writer.addData4(1, 0, 0, 1)  # Red
        
        lines = GeomLines(Geom.UHDynamic)
        for i in range(len(path) - 1):
            lines.addVertices(i, i + 1)
        
        geom = Geom(vdata)
        geom.addPrimitive(lines)
        
        node = GeomNode('ray_lines')
        node.addGeom(geom)
        
        self.ray_line_node = self.render.attachNewNode(node)
        self.ray_line_node.setRenderModeThickness(3)
        
    def check_wire_hit(self, ray_pos, ray_dir):
        """Check if ray enters the wire entrance"""
        position = self.wire.getPythonTag('position')
        direction = self.wire.getPythonTag('direction')
        radius = self.wire.getPythonTag('radius')
        
        # Update position
        position = self.wire.getPos()
        
        to_wire = position - ray_pos
        dist = to_wire.length()
        
        dot_product = ray_dir.dot(to_wire)
        if dot_product <= 0 or dist > 200:
            return None
        
        wire_dot = ray_dir.dot(direction)
        if abs(wire_dot) < 0.1:
            return None
        
        plane_dist = to_wire.dot(direction)
        t = plane_dist / wire_dot
        if t < 0:
            return None
        
        intersection = ray_pos + ray_dir * t
        to_intersection = intersection - position
        axial_dist = to_intersection.dot(direction)
        radial = to_intersection - direction * axial_dist
        radial_dist = radial.length()
        
        if radial_dist <= radius:
            return intersection, radial_dist
        return None
    
    def calculate_intensity(self, hit_info, ray_dir):
        """Calculate coupling intensity"""
        if not hit_info:
            return 0.0  # Minimal intensity for miss
        
        intersection, radial_dist = hit_info
        radius = self.wire.getPythonTag('radius')
        direction = self.wire.getPythonTag('direction')
        
        eta_lateral = math.exp(-(radial_dist ** 2) / (radius ** 2))
        alignment = abs(ray_dir.dot(direction))
        eta_angular = alignment ** 2

        noise = np.random.normal(-.01, 0.002)
        
        return eta_lateral * eta_angular #+ noise
    
    def evaluate_intensity_for_angles(self, angles):
        """Evaluate intensity for given mirror angles [ax1, ay1, ax2, ay2]

        This is the objective function that optimizers can call.
        It sets the disk orientations, traces the ray, and returns intensity.
        The disks remain at the evaluated position (for visual feedback).
        """
        ax1, ay1, ax2, ay2 = angles

        # Create rotation quaternions for disk 1
        rot1_x = LQuaternionf()
        rot1_x.setFromAxisAngle(math.degrees(ax1), LVector3f(1, 0, 0))
        rot1_y = LQuaternionf()
        rot1_y.setFromAxisAngle(math.degrees(ay1), LVector3f(0, 1, 0))
        
        # Create rotation quaternions for disk 2
        rot2_x = LQuaternionf()
        rot2_x.setFromAxisAngle(math.degrees(ax2), LVector3f(1, 0, 0))
        rot2_y = LQuaternionf()
        rot2_y.setFromAxisAngle(math.degrees(ay2), LVector3f(0, 1, 0))
        
        # Get initial orientations from optimizer
        if hasattr(self.optimizer, 'initial_quat1'):
            initial_quat1 = self.optimizer.initial_quat1
            initial_quat2 = self.optimizer.initial_quat2
        else:
            initial_quat1 = self.disk1_initial_quat
            initial_quat2 = self.disk2_initial_quat
        
        # Combine with initial orientations
        quat1 = rot1_y * rot1_x * initial_quat1
        quat2 = rot2_y * rot2_x * initial_quat2
        
        # Set disk orientations (keep them visible at this position)
        self.disk1.setQuat(quat1)
        self.disk2.setQuat(quat2)
        
        # Update normal vectors
        forward = LVector3f(0, 0, 1)
        new_normal1 = self.disk1.getQuat().xform(forward)
        new_normal1.normalize()
        self.disk1.setPythonTag('normal', new_normal1)
        
        new_normal2 = self.disk2.getQuat().xform(forward)
        new_normal2.normalize()
        self.disk2.setPythonTag('normal', new_normal2)
        
        # Trace ray and calculate intensity
        path, final_pos, final_dir = self.trace_ray_3d()
        hit_info = self.check_wire_hit(final_pos, final_dir)
        intensity = self.calculate_intensity(hit_info, final_dir)
        
        # Redraw the ray path
        self.draw_ray_path(path)
        # Update intensity text during optimization
        self.intensity_text.setText(f'Wire Coupling Intensity: {intensity:.4f}')
        
        # Notify optimizer about the actual angles used for this evaluation
        if hasattr(self, 'optimizer') and self.optimizer:
            self.optimizer.last_angles = angles
            if hasattr(self.optimizer, 'best_intensity') and intensity > self.optimizer.best_intensity:
                self.optimizer.best_intensity = intensity
        
        return intensity
    
    def update_task(self, task):
        """Main update loop"""
        angle_step = 1.5 / 25
        
        # Handle disk 1 rotation (arrow keys)
        if self.keys.get('left'):
            self.rotate_disk(self.disk1, 'z', angle_step)
        if self.keys.get('right'):
            self.rotate_disk(self.disk1, 'z', -angle_step)
        if self.keys.get('up'):
            self.rotate_disk(self.disk1, 'x', angle_step)
        if self.keys.get('down'):
            self.rotate_disk(self.disk1, 'x', -angle_step)
        
        # Handle disk 2 rotation (WASD)
        if self.keys.get('a'):
            self.rotate_disk(self.disk2, 'z', angle_step)
        if self.keys.get('d'):
            self.rotate_disk(self.disk2, 'z', -angle_step)
        if self.keys.get('w'):
            self.rotate_disk(self.disk2, 'y', angle_step)
        if self.keys.get('s'):
            self.rotate_disk(self.disk2, 'y', -angle_step)

        # Handle camera orbit with mouse
        if self.mouseWatcherNode.hasMouse():
            mouse_x = self.mouseWatcherNode.getMouseX()
            mouse_y = self.mouseWatcherNode.getMouseY()
            
            if self.is_orbiting:
                dx = mouse_x - self.last_mouse_x
                dy = mouse_y - self.last_mouse_y
                
                if abs(dx) > 0.001 or abs(dy) > 0.001:
                    # Rotate camera around target
                    self.camera_heading -= dx * 100  # Horizontal rotation
                    self.camera_pitch += dy * 100    # Vertical rotation
                    
                    # Clamp pitch to avoid gimbal lock
                    self.camera_pitch = max(-89, min(89, self.camera_pitch))
                    
                    self.update_camera_position()
                
                self.last_mouse_x = mouse_x
                self.last_mouse_y = mouse_y
                
            elif self.is_panning:
                dx = mouse_x - self.last_mouse_x
                dy = mouse_y - self.last_mouse_y
                
                if abs(dx) > 0.001 or abs(dy) > 0.001:
                    # Pan camera target
                    # Get camera's right and up vectors
                    heading_rad = math.radians(self.camera_heading)
                    
                    # Right vector (perpendicular to heading in XY plane)
                    right_x = math.cos(heading_rad)
                    right_y = -math.sin(heading_rad)
                    
                    # Up vector (always Z in world space for simpler panning)
                    pan_speed = self.camera_distance * 0.5
                    
                    self.camera_target.x -= right_x * dx * pan_speed
                    self.camera_target.y -= right_y * dx * pan_speed
                    self.camera_target.z += dy * pan_speed
                    
                    self.update_camera_position()
                
                self.last_mouse_x = mouse_x
                self.last_mouse_y = mouse_y
        
        # Handle optimizer
        if self.optimizer.is_running:
            # Special handling for NelderMeadOptimizer (runs scipy synchronously)
            if isinstance(self.optimizer, ScipyOptimizer):
                # NelderMead has already run scipy.optimize.minimize in start()
                # Just display the best result found
                if self.optimizer.is_complete():
                    best_state = self.optimizer.get_best_state()
                    if best_state:
                        quat1, quat2 = best_state
                        self.disk1.setQuat(quat1)
                        self.disk2.setQuat(quat2)
                        
                        # Update normal vectors
                        forward = LVector3f(0, 0, 1)
                        new_normal1 = self.disk1.getQuat().xform(forward)
                        new_normal1.normalize()
                        self.disk1.setPythonTag('normal', new_normal1)
                        
                        new_normal2 = self.disk2.getQuat().xform(forward)
                        new_normal2.normalize()
                        self.disk2.setPythonTag('normal', new_normal2)
                    
                    # Update display
                    samples = len(self.optimizer.history)
                    opt_name = "Nelder-Mead"
                    self.optimizer_text.setText(f'{opt_name}: COMPLETE | Samples: {samples} | Best: {self.optimizer.best_intensity:.4f}')
                else:
                    # Still running (shouldn't happen with synchronous scipy, but keep for safety)
                    samples = len(self.optimizer.history)
                    opt_name = "Nelder-Mead"
                    self.optimizer_text.setText(f'{opt_name}: RUNNING | Samples: {samples} | Best: {self.optimizer.best_intensity:.4f}')
            else:
                # Other optimizers use async get_next_state() pattern
                state = self.optimizer.get_next_state()
                if state:
                    quat1, quat2 = state
                    self.disk1.setQuat(quat1)
                    self.disk2.setQuat(quat2)
                    
                    # Update normal vectors
                    forward = LVector3f(0, 0, 1)
                    new_normal1 = self.disk1.getQuat().xform(forward)
                    new_normal1.normalize()
                    self.disk1.setPythonTag('normal', new_normal1)
                    
                    new_normal2 = self.disk2.getQuat().xform(forward)
                    new_normal2.normalize()
                    self.disk2.setPythonTag('normal', new_normal2)
                
                # Check if search is complete or show progress
                if isinstance(self.optimizer, GaussianProcessOptimizer):
                    samples = self.optimizer.current_sample
                    phase = "Explore" if samples < self.optimizer.n_initial_samples else "Exploit"
                    status = "COMPLETE" if self.optimizer.is_complete() else "RUNNING"
                    self.optimizer_text.setText(f'GP: {status} | {phase} | {samples}/{self.optimizer.n_total_samples} | Best: {self.optimizer.best_intensity:.4f}')
                elif self.optimizer.is_complete():
                    best_state = self.optimizer.get_best_state()
                    if best_state:
                        quat1, quat2 = best_state
                        self.disk1.setQuat(quat1)
                        self.disk2.setQuat(quat2)
                        
                        # Update normal vectors
                        forward = LVector3f(0, 0, 1)
                        new_normal1 = self.disk1.getQuat().xform(forward)
                        new_normal1.normalize()
                        self.disk1.setPythonTag('normal', new_normal1)
                        
                        new_normal2 = self.disk2.getQuat().xform(forward)
                        new_normal2.normalize()
                        self.disk2.setPythonTag('normal', new_normal2)
                    # For IterativeGridSearch, regenerate grid
                    if hasattr(self.optimizer, 'generate_search_grid'):
                        self.optimizer.generate_search_grid(self.disk1, self.disk2)
                    progress = min(100, (self.optimizer.best_intensity * 100))
                    
                    if hasattr(self.optimizer, 'samples_taken'):
                        samples = self.optimizer.samples_taken
                        self.optimizer_text.setText(f'Grid Search: Samples: {samples} | Best: {self.optimizer.best_intensity:.4f} ({progress:.0f}%)')
        
        # Trace and visualize ray
        path, final_pos, final_dir = self.trace_ray_3d()
        self.draw_ray_path(path)
        
        # Calculate intensity
        hit_info = self.check_wire_hit(final_pos, final_dir)
        intensity = self.calculate_intensity(hit_info, final_dir)
        
        self.intensity_text.setText(f'Wire Coupling Intensity: {intensity:.4f}')
        
        # Update optimizer with current intensity (for async optimizers only)
        if self.optimizer.is_running and not isinstance(self.optimizer, ScipyOptimizer):
            if hasattr(self.optimizer, 'update'):
                self.optimizer.update(intensity, self.disk1, self.disk2)
        
        return Task.cont

# Run the simulation
if __name__ == '__main__':
    app = LaserSimulation()
    app.run()
