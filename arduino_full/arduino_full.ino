/*
 * Arduino Nelder–Mead 3D Laser Alignment Optimizer
 * Controls four 2048-step steppers (e.g., 28BYJ-48) for two 2-axis mirrors.
 * Uses 4D Nelder-Mead optimization to find mirror angles that maximize 
 * laser coupling intensity into a fiber optic or photodetector.
 * 
 * Hardware:
 * - 4x Stepper motors (28BYJ-48 with ULN2003 drivers)
 * - 1x Photodetector (analog output on A0)
 * - 2x Gimbal-mounted mirrors (2 steppers each for X/Y control)
 * 
 * Optimization:
 * - Algorithm: Nelder-Mead simplex method (derivative-free)
 * - Dimensions: 4D (ax1, ay1, ax2, ay2)
 * - Typical convergence: 50-200 iterations
 */

#include <Arduino.h>
#include <Stepper.h>

/* ------------------- CONFIGURATION ------------------- */

// Stepper configuration
#define STEPS_PER_REV 2048
#define REV_PER_DEGREE 2   // Product page lists 0.5º per revolution or 2 rev/degree

// Stepper pin mapping TODO: need to adjust this to fit the physical board
Stepper stepper_m1x(STEPS_PER_REV, 8, 10, 9, 11);
Stepper stepper_m1y(STEPS_PER_REV, 4, 6, 5, 7);
Stepper stepper_m2x(STEPS_PER_REV, 22, 24, 23, 25);
Stepper stepper_m2y(STEPS_PER_REV, 18, 20, 19, 21);

static long stepper_m1x_current_pos = 0;
static long stepper_m1y_current_pos = 0;
static long stepper_m2x_current_pos = 0;
static long stepper_m2y_current_pos = 0;

// Sensor input
#define SENSOR_PIN A0
#define ADC_SAMPLES 8
#define SETTLE_DELAY_MS 50

// Optimization parameters (in steps)
#define INITIAL_STEP_SIZE_STEPS 100        // initial simplex step size (motor steps)
#define MAX_ITERATIONS 500
#define REFLECTION_COEFFICIENT 1.0         // alpha: reflection
#define EXPANSION_COEFFICIENT 2.0          // gamma: expansion
#define CONTRACTION_COEFFICIENT 0.5        // rho: contraction
#define SHRINK_COEFFICIENT 0.5             // sigma: shrink
#define FUNCTION_TOLERANCE 1e-8            // Very tight convergence criterion
#define PARAMETER_TOLERANCE 1.0            // tolerance in steps

// Testing mode
bool SIMULATED = false;  // Set to true to use simulated intensity (no motor movement)

/* ------------------- MIRROR / STEPPER UTILS ------------------- */

/**
 * Move a single stepper motor to a target position in steps.
 * 
 * Calculates the relative movement needed from the current position,
 * executes the move, and updates the position tracker. 
 * Note: stepper.step() is BLOCKING and will halt execution until complete.
 * 
 * @param stepper Reference to the Stepper object to control
 * @param current_position Reference to the position tracker (will be updated)
 * @param target_steps Target position in motor steps
 */
void setMirrorPosition(Stepper &stepper, long &current_position, long target_steps) {
  long step_difference = target_steps - current_position;
  stepper.step(step_difference);
  current_position = target_steps;
}

/**
 * Move all four stepper motors to their target positions.
 * 
 * Controls the X and Y axes for both mirrors (M1 and M2) by moving each
 * stepper motor sequentially. Motors move one at a time (not simultaneously)
 * due to the blocking nature of the Stepper library.
 * 
 * @param mirror1_x_steps Mirror 1 X-axis position in steps
 * @param mirror1_y_steps Mirror 1 Y-axis position in steps
 * @param mirror2_x_steps Mirror 2 X-axis position in steps
 * @param mirror2_y_steps Mirror 2 Y-axis position in steps
 */
void setMirrorPositions(long mirror1_x_steps, long mirror1_y_steps, 
                        long mirror2_x_steps, long mirror2_y_steps) {
  if (SIMULATED) {
    // In simulation mode, just update position trackers without moving motors
    stepper_m1x_current_pos = mirror1_x_steps;
    stepper_m1y_current_pos = mirror1_y_steps;
    stepper_m2x_current_pos = mirror2_x_steps;
    stepper_m2y_current_pos = mirror2_y_steps;
    return;
  }
  
  setMirrorPosition(stepper_m1x, stepper_m1x_current_pos, mirror1_x_steps);
  setMirrorPosition(stepper_m1y, stepper_m1y_current_pos, mirror1_y_steps);
  setMirrorPosition(stepper_m2x, stepper_m2x_current_pos, mirror2_x_steps);
  setMirrorPosition(stepper_m2y, stepper_m2y_current_pos, mirror2_y_steps);
}

/* ------------------- MEASUREMENT ------------------- */

/**
 * Measure light intensity at a given mirror configuration.
 * 
 * Moves all four mirrors to the specified positions, waits for mechanical settling,
 * then takes multiple ADC readings from the photodetector and returns their
 * average. Higher values indicate brighter light (better alignment).
 * 
 * @param mirror1_x_steps Mirror 1 X-axis position in steps
 * @param mirror1_y_steps Mirror 1 Y-axis position in steps
 * @param mirror2_x_steps Mirror 2 X-axis position in steps
 * @param mirror2_y_steps Mirror 2 Y-axis position in steps
 * @return Average ADC reading (0-1023 for 10-bit ADC), normalized to [0, 1]
 */
double measureIntensity(long mirror1_x_steps, long mirror1_y_steps, 
                        long mirror2_x_steps, long mirror2_y_steps) {
  if (SIMULATED) {
    return simulatedMeasureIntensity(mirror1_x_steps, mirror1_y_steps, 
                                      mirror2_x_steps, mirror2_y_steps);
  }
  
  setMirrorPositions(mirror1_x_steps, mirror1_y_steps, mirror2_x_steps, mirror2_y_steps);
  delay(SETTLE_DELAY_MS);

  long sum = 0;
  for (int i = 0; i < ADC_SAMPLES; i++) {
    sum += analogRead(SENSOR_PIN);
  }
  double average = (double)sum / ADC_SAMPLES;
  return average / 1023.0; // Normalize to [0, 1]
}

double simulatedMeasureIntensity(long m1x, long m1y, long m2x, long m2y) {
  // Simulated "optimal" positions (in steps)
  const long OPTIMAL_M1X = 500;
  const long OPTIMAL_M1Y = 300;
  const long OPTIMAL_M2X = 700;
  const long OPTIMAL_M2Y = 400;

  // Calculate squared distance from optimal
  long dx = m1x - OPTIMAL_M1X;
  long dy = m1y - OPTIMAL_M1Y;
  long dz = m2x - OPTIMAL_M2X;
  long dw = m2y - OPTIMAL_M2Y;
  long distance_squared = dx*dx + dy*dy + dz*dz + dw*dw;

  // Simulated intensity: max at optimal, decreases with distance
  double intensity = exp(-0.0001 * distance_squared); // Gaussian falloff
  return intensity;
}

/* ------------------- SIMPLE UTILS ------------------- */

/**
 * Copy values from one vector to another.
 * 
 * Simple utility for vector operations in the Nelder-Mead algorithm.
 * 
 * @param destination Destination array pointer
 * @param source Source array pointer
 * @param num_elements Number of elements to copy
 */
void copyVector(double *destination, const double *source, int num_elements) {
  for (int i = 0; i < num_elements; i++) {
    destination[i] = source[i];
  }
}

/**
 * Calculate the centroid of n points in N-dimensional space.
 * 
 * Used in the Nelder-Mead algorithm to find the center of the best n
 * vertices of the simplex before performing reflection/expansion operations.
 * 
 * @param centroid Output: centroid coordinates
 * @param points Array of n points, each with dimensions
 * @param num_points Number of points
 * @param dimensions Dimensionality of each point
 */
void calculateCentroid(double *centroid, double **points, int num_points, int dimensions) {
  for (int dimension = 0; dimension < dimensions; dimension++) {
    centroid[dimension] = 0.0;
    for (int i = 0; i < num_points; i++) {
      centroid[dimension] += points[i][dimension];
    }
    centroid[dimension] /= (double)num_points;
  }
}

/* ------------------- NELDER–MEAD (4D) ------------------- */

/**
 * Nelder-Mead optimization algorithm for 4D parameter space (in motor steps).
 * 
 * A derivative-free optimization method that uses a simplex (5 points in 4D)
 * to iteratively search for the motor positions (in steps) that maximize light intensity.
 * The algorithm reflects, expands, contracts, or shrinks the simplex based on
 * function evaluations at each vertex.
 * 
 * The algorithm minimizes -intensity (equivalent to maximizing intensity).
 * Operations performed at each iteration:
 * - Reflection: Try a point mirrored across the centroid from the worst point
 * - Expansion: If reflection is good, try going further in that direction
 * - Contraction: If reflection is bad, try a point between centroid and worst
 * - Shrink: If contraction fails, pull all points toward the best point
 * 
 * @param initial_guess Initial guess for parameters [m1x, m1y, m2x, m2y] in steps
 * @param step_size Initial simplex size in steps
 * @param max_iterations Maximum number of iterations before termination
 * @param function_tolerance Function value tolerance for convergence (checks std dev of function values)
 * @param parameter_tolerance Parameter tolerance for convergence (currently unused)
 * @param best_parameters_out Output: best parameters found [m1x, m1y, m2x, m2y] in steps
 * @param best_intensity_out Output: best intensity value achieved
 */
void nelderMead4D(long initial_guess[4], long step_size, int max_iterations,
                  double function_tolerance, double parameter_tolerance,
                  long best_parameters_out[4], double *best_intensity_out) {

  const int NUM_DIMENSIONS = 4;
  const int NUM_VERTICES = NUM_DIMENSIONS + 1;  // 5 vertices for 4D simplex

  // Simplex vertices (5 points in 4D) - using doubles for intermediate calculations
  double simplex_vertices[NUM_VERTICES][NUM_DIMENSIONS];              // The coordinates of each vertex
  double function_values[NUM_VERTICES];                               // The function values at each vertex

  // Build initial simplex using regular simplex construction
  // Vertex 0: initial guess
  for (int dimension = 0; dimension < NUM_DIMENSIONS; dimension++) {
    simplex_vertices[0][dimension] = (double)initial_guess[dimension];
  }
  
  // Vertices 1-4: perturb each dimension
  for (int vertex = 1; vertex <= NUM_DIMENSIONS; vertex++) {
    for (int dimension = 0; dimension < NUM_DIMENSIONS; dimension++) {
      simplex_vertices[vertex][dimension] = (double)initial_guess[dimension];
    }
    simplex_vertices[vertex][vertex-1] += (double)step_size;  // Perturb dimension vertex-1
  }

  // Measure initial points
  for (int vertex = 0; vertex < NUM_VERTICES; vertex++) {
    function_values[vertex] = -measureIntensity(
      lround(simplex_vertices[vertex][0]), 
      lround(simplex_vertices[vertex][1]), 
      lround(simplex_vertices[vertex][2]), 
      lround(simplex_vertices[vertex][3])
    ); // negative => minimize
  }

  int iteration = 0;
  while (iteration < max_iterations) {
    // Sort vertices by function_values ascending (best first)
    for (int i = 0; i < NUM_VERTICES - 1; i++) {
      for (int j = i + 1; j < NUM_VERTICES; j++) {
        if (function_values[j] < function_values[i]) {
          double temp_value = function_values[i]; 
          function_values[i] = function_values[j]; 
          function_values[j] = temp_value;
          
          double temp_vertex[NUM_DIMENSIONS]; 
          copyVector(temp_vertex, simplex_vertices[i], NUM_DIMENSIONS);
          copyVector(simplex_vertices[i], simplex_vertices[j], NUM_DIMENSIONS); 
          copyVector(simplex_vertices[j], temp_vertex, NUM_DIMENSIONS);
        }
      }
    }

    // Compute centroid of best n points (exclude worst)
    double centroid[NUM_DIMENSIONS];
    double *best_points[NUM_DIMENSIONS];
    for (int i = 0; i < NUM_DIMENSIONS; i++) {
      best_points[i] = simplex_vertices[i];
    }
    calculateCentroid(centroid, best_points, NUM_DIMENSIONS, NUM_DIMENSIONS);

    // Reflection
    double reflected_point[NUM_DIMENSIONS];
    for (int dimension = 0; dimension < NUM_DIMENSIONS; dimension++) {
      reflected_point[dimension] = centroid[dimension] + 
        REFLECTION_COEFFICIENT * (centroid[dimension] - simplex_vertices[NUM_VERTICES-1][dimension]);
    }
    double reflected_value = -measureIntensity(
      lround(reflected_point[0]), lround(reflected_point[1]), 
      lround(reflected_point[2]), lround(reflected_point[3])
    );

    if (reflected_value < function_values[0]) {
      // Expansion
      double expanded_point[NUM_DIMENSIONS];
      for (int dimension = 0; dimension < NUM_DIMENSIONS; dimension++) {
        expanded_point[dimension] = centroid[dimension] + 
          EXPANSION_COEFFICIENT * (reflected_point[dimension] - centroid[dimension]);
      }
      double expanded_value = -measureIntensity(
        lround(expanded_point[0]), lround(expanded_point[1]), 
        lround(expanded_point[2]), lround(expanded_point[3])
      );
      if (expanded_value < reflected_value) {
        copyVector(simplex_vertices[NUM_VERTICES-1], expanded_point, NUM_DIMENSIONS);
        function_values[NUM_VERTICES-1] = expanded_value;
      } else {
        copyVector(simplex_vertices[NUM_VERTICES-1], reflected_point, NUM_DIMENSIONS);
        function_values[NUM_VERTICES-1] = reflected_value;
      }
    } else if (reflected_value < function_values[NUM_DIMENSIONS-1]) {
      copyVector(simplex_vertices[NUM_VERTICES-1], reflected_point, NUM_DIMENSIONS);
      function_values[NUM_VERTICES-1] = reflected_value;
    } else {
      // Contraction
      double contracted_point[NUM_DIMENSIONS];
      if (reflected_value < function_values[NUM_VERTICES-1]) {
        // Outside contraction
        for (int dimension = 0; dimension < NUM_DIMENSIONS; dimension++) {
          contracted_point[dimension] = centroid[dimension] + 
            CONTRACTION_COEFFICIENT * (reflected_point[dimension] - centroid[dimension]);
        }
      } else {
        // Inside contraction
        for (int dimension = 0; dimension < NUM_DIMENSIONS; dimension++) {
          contracted_point[dimension] = centroid[dimension] + 
            CONTRACTION_COEFFICIENT * (simplex_vertices[NUM_VERTICES-1][dimension] - centroid[dimension]);
        }
      }
      double contracted_value = -measureIntensity(
        lround(contracted_point[0]), lround(contracted_point[1]), 
        lround(contracted_point[2]), lround(contracted_point[3])
      );
      if (contracted_value < function_values[NUM_VERTICES-1]) {
        copyVector(simplex_vertices[NUM_VERTICES-1], contracted_point, NUM_DIMENSIONS);
        function_values[NUM_VERTICES-1] = contracted_value;
      } else {
        // Shrink towards best
        for (int vertex = 1; vertex < NUM_VERTICES; vertex++) {
          for (int dimension = 0; dimension < NUM_DIMENSIONS; dimension++) {
            simplex_vertices[vertex][dimension] = simplex_vertices[0][dimension] + 
              SHRINK_COEFFICIENT * (simplex_vertices[vertex][dimension] - simplex_vertices[0][dimension]);
          }
          function_values[vertex] = -measureIntensity(
            lround(simplex_vertices[vertex][0]), lround(simplex_vertices[vertex][1]), 
            lround(simplex_vertices[vertex][2]), lround(simplex_vertices[vertex][3])
          );
        }
      }
    }

    // Convergence check (std dev of function values)
    double mean_value = 0.0;
    for (int i = 0; i < NUM_VERTICES; i++) {
      mean_value += function_values[i];
    }
    mean_value /= NUM_VERTICES;
    
    double variance = 0.0;
    for (int i = 0; i < NUM_VERTICES; i++) {
      double difference = function_values[i] - mean_value;
      variance += difference * difference;
    }
    variance /= NUM_VERTICES;
    double standard_deviation = sqrt(variance);
    
    // Progress report every 10 iterations (or first few iterations for debugging)
    if (iteration % 10 == 0) {
      Serial.print("Iteration "); Serial.print(iteration);
      Serial.print(": best_intensity = "); Serial.println(-function_values[0], 4);
    }
    
    // Only check convergence after some iterations and when we have meaningful values
    if (iteration > 10 && standard_deviation < function_tolerance && mean_value < -0.01) {
      Serial.print("Converged at iteration ");
      Serial.println(iteration);
      break;
    }

    iteration++;
  }

  // Convert best parameters from double to long (round to nearest step)
  for (int dimension = 0; dimension < NUM_DIMENSIONS; dimension++) {
    best_parameters_out[dimension] = lround(simplex_vertices[0][dimension]);
  }
  *best_intensity_out = -function_values[0]; // convert back to positive intensity
}

/* ------------------- ARDUINO MAIN ------------------- */

void setup() {
  Serial.begin(115200);
  analogReference(DEFAULT);
  Serial.println("Starting 4D mirror alignment optimizer...");
  Serial.println("Controlling 4 stepper motors for 3D laser alignment");
  Serial.println("Optimization units: motor steps");

  // Check if we're running in simulation mode
  if (SIMULATED) {
    Serial.println("\n*** SIMULATION MODE ENABLED ***");
    Serial.println("Motors will NOT move. Using simulated intensity function.");
    Serial.println("Optimal position: M1X=500, M1Y=300, M2X=700, M2Y=400\n");
  }

  // Configure stepper speeds (RPM)
  stepper_m1x.setSpeed(10);
  stepper_m1y.setSpeed(10);
  stepper_m2x.setSpeed(10);
  stepper_m2y.setSpeed(10);

  // Initial guess: all positions at 0 steps (home position)
  long initial_guess[4] = {0, 0, 0, 0};  // [m1x, m1y, m2x, m2y] in steps
  long best_parameters[4];
  double best_intensity;

  Serial.println("Starting Nelder-Mead optimization...");
  Serial.print("Initial simplex step size: "); 
  Serial.print(INITIAL_STEP_SIZE_STEPS); 
  Serial.println(" steps");
  
  unsigned long start_time = millis();
  nelderMead4D(initial_guess, INITIAL_STEP_SIZE_STEPS, MAX_ITERATIONS, 
               FUNCTION_TOLERANCE, PARAMETER_TOLERANCE, best_parameters, &best_intensity);
  unsigned long elapsed_time = millis() - start_time;

  Serial.println("\n=== Optimization Complete ===");
  Serial.print("Elapsed time: "); 
  Serial.print(elapsed_time / 1000.0, 2); 
  Serial.println(" seconds");
  Serial.print("Best intensity = "); Serial.println(best_intensity, 4);
  Serial.println("Best positions (steps):");
  Serial.print("  Mirror 1 X: "); Serial.println(best_parameters[0]);
  Serial.print("  Mirror 1 Y: "); Serial.println(best_parameters[1]);
  Serial.print("  Mirror 2 X: "); Serial.println(best_parameters[2]);
  Serial.print("  Mirror 2 Y: "); Serial.println(best_parameters[3]);
  
  Serial.println("Best positions (degrees):");
  Serial.print("  Mirror 1 X: "); 
  Serial.println((double)best_parameters[0] / (STEPS_PER_REV * REV_PER_DEGREE), 3);
  Serial.print("  Mirror 1 Y: "); 
  Serial.println((double)best_parameters[1] / (STEPS_PER_REV * REV_PER_DEGREE), 3);
  Serial.print("  Mirror 2 X: "); 
  Serial.println((double)best_parameters[2] / (STEPS_PER_REV * REV_PER_DEGREE), 3);
  Serial.print("  Mirror 2 Y: "); 
  Serial.println((double)best_parameters[3] / (STEPS_PER_REV * REV_PER_DEGREE), 3);

  if (SIMULATED) {
    // In simulation mode, verify we found the correct optimal position
    Serial.println("\n=== Simulation Test Results ===");
    Serial.println("Expected optimal position: M1X=500, M1Y=300, M2X=700, M2Y=400");
    Serial.print("Distance from optimal: ");
    long dx = best_parameters[0] - 500;
    long dy = best_parameters[1] - 300;
    long dz = best_parameters[2] - 700;
    long dw = best_parameters[3] - 400;
    double distance = sqrt(dx*dx + dy*dy + dz*dz + dw*dw);
    Serial.print(distance, 2);
    Serial.println(" steps");
    
    if (distance < 10.0) {
      Serial.println("TEST PASSED: Found optimal position within 10 steps!");
    } else {
      Serial.println("TEST WARNING: Optimal position not found accurately.");
    }
  } else {
    // Move to best configuration
    Serial.println("\nMoving to optimal configuration...");
    setMirrorPositions(best_parameters[0], best_parameters[1], best_parameters[2], best_parameters[3]);
    Serial.println("Done!");
  }
}

void loop() {
  // nothing
}
