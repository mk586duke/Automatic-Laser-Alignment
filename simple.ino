#include "Stepper.h"

const int stepsPerRevolution = 2048; // value from datasheet

// Wiring:
// Pin 22 to IN1 on the ULN2003 driver
// Pin 24 to IN2 on the ULN2003 driver
// Pin 26 to IN3 on the ULN2003 driver
// Pin 28 to IN4 on the ULN2003 driver

// Create stepper object called 'myStepper', note the pin order:
// Stepper myStepper = Stepper(stepsPerRevolution, 22, 26, 24, 28);
Stepper stepper_m1x(STEPS_PER_REV, 8, 10, 9, 11);
Stepper stepper_m1y(STEPS_PER_REV, 4, 6, 5, 7);
Stepper stepper_m2x(STEPS_PER_REV, 22, 24, 23, 25);
Stepper stepper_m2y(STEPS_PER_REV, 18, 20, 19, 21);

// Simple test program to rotate the stepper motor into position
Stepper activeStepper = stepper_m1x;
int rpm = 10;
int direction = 1; // 1 for clockwise, -1 for counterclockwise

void setup() {
  activeStepper.setSpeed(rpm);
  // Tested maximum 17 rpm
  // Function minimum 1 rpm, can reduce much more with delays in loop
  // Begin Serial communication at a baud rate of 9600:
  Serial.begin(9600);
}

void loop() {
  // Step one revolution 2048 steps:
  Serial.println("Stepping 1 revolution...");
  activeStepper.step(stepsPerRevolution * direction);
  delay(500);
}
