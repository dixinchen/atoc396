int value = 100;

void setup() {
  // put your setup code here, to run once:
//  TCCR2B = TCCR2B & 0b11111000 | 0x06;
//  TCCR1B = TCCR1B & 0b11111000 | 0x04;
  pinMode(3, OUTPUT);
  analogWrite(3, value);
  pinMode(5, OUTPUT);
  analogWrite(5, value);
  pinMode(9, OUTPUT);
  analogWrite(9, value);
  pinMode(10, OUTPUT);
  analogWrite(10, value);
}

void loop() {
}
