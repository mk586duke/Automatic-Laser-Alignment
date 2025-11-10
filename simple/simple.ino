#include "Stepper.h"

const int stepsPerRevolution = 2048; // value from datasheet

// Wiring:
// m1x: in1=14, in2=15, in3=16, in4=17
// m1y: in1=21, in2=20, in3=19, in4=18
// m2x: in1=22, in2=24, in3=23, in4=25
// m2y: in1=0, in2=1, in3=2, in4=3

// Create stepper object called 'myStepper', note the pin order:
// Stepper myStepper = Stepper(stepsPerRevolution, 22, 26, 24, 28);
int m1y_in1 = 14; int m1y_in2 = 15; int m1y_in3 = 16; int m1y_in4 = 17;
int m1x_in1 = 21; int m1x_in2 = 20; int m1x_in3 = 19; int m1x_in4 = 18;
int m2y_in1 = 7; int m2y_in2 = 6; int m2y_in3 = 5; int m2y_in4 = 4;
int m2x_in1 = 10; int m2x_in2 = 11; int m2x_in3 = 12; int m2x_in4 = 13;
Stepper stepper_m1x(stepsPerRevolution, m1x_in1, m1x_in3, m1x_in2, m1x_in4);
Stepper stepper_m1y(stepsPerRevolution, m1y_in1, m1y_in3, m1y_in2, m1y_in4);
Stepper stepper_m2x(stepsPerRevolution, m2x_in1, m2x_in3, m2x_in2, m2x_in4);
Stepper stepper_m2y(stepsPerRevolution, m2y_in1, m2y_in3, m2y_in2, m2y_in4);

// Simple test program to rotate the stepper motor into position
Stepper activeStepper = stepper_m2y;
int rpm = 10;
int direction = 0; // 1 for clockwise, -1 for counterclockwise

void setup() {
  activeStepper.setSpeed(rpm);
  // Tested maximum 17 rpm
  // Function minimum 1 rpm, can reduce much more with delays in loop
  // Begin Serial communication at a baud rate of 9600:
  Serial.begin(9600);
  delay(1000);
  Serial.println("Arduino started!");
  Serial.println("Motor: m1y, RPM: 5, Direction: -1");
}

void loop() {
  // stepper_m1x.setSpeed(rpm);
  Serial.println("Stepping 1 revolution!...");
  activeStepper.step(stepsPerRevolution * direction);
  // Motors stopped
  delay(1000);
}
