#include "Stepper.h"

const int stepsPerRevolution = 2048; // value from datasheet

// Wiring:
// Pin 22 to IN1 on the ULN2003 driver
// Pin 24 to IN2 on the ULN2003 driver
// Pin 26 to IN3 on the ULN2003 driver
// Pin 28 to IN4 on the ULN2003 driver

// Create stepper object called 'myStepper', note the pin order:
Stepper myStepper = Stepper(stepsPerRevolution, 22, 26, 24, 28);

void setup() {
  myStepper.setSpeed(4);
  // Tested maximum 17 rpm
  // Function minimum 1 rpm, can reduce much more with delays in loop
  // Begin Serial communication at a baud rate of 9600:
  Serial.begin(9600);
}

void loop() {
  // Step one revolution in one direction:
  Serial.println("clockwise");
  myStepper.step(stepsPerRevolution);
  delay(500);
  
  // Step one revolution in the other direction:
  Serial.println("counterclockwise");
  myStepper.step(-stepsPerRevolution);
  delay(500);
}
