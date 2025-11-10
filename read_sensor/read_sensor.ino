/*
 * Simple Analog Sensor Reader
 * Reads from pin A0 and prints the value to Serial
 */

#define SENSOR_PIN A0
#define BAUD_RATE 9600
#define SAMPLE_DELAY_MS 100  // Delay between readings (100ms = 10 readings/sec)

void setup() {
  Serial.begin(BAUD_RATE);
  analogReference(DEFAULT);
  
  delay(1000);
  Serial.println("Analog Sensor Reader");
  Serial.println("Reading from pin A0");
  Serial.println("Value range: 0-1023 (10-bit ADC)");
  Serial.println("---");
}

void loop() {
  int raw_value = analogRead(SENSOR_PIN);
  double voltage = (raw_value / 1023.0) * 5.0;  // Assuming 5V reference
  
  Serial.print("Raw: ");
  Serial.print(raw_value);
  Serial.print(" | Voltage: ");
  Serial.print(voltage, 3);
  Serial.println(" V");
  
  delay(SAMPLE_DELAY_MS);
}
