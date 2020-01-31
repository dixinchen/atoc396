#define BLUE 3
#define GREEN 5
#define RED 6

void setup()
{
pinMode(RED, OUTPUT);
pinMode(GREEN, OUTPUT);
pinMode(BLUE, OUTPUT);
Serial.begin(9600);
}

// define variables
int redValue;
int greenValue;
int blueValue;

// main loop
void loop()
{
//analogWrite(RED, 240);
//analogWrite(GREEN, 30);
//analogWrite(BLUE, 10);
//delay(1000);
//analogWrite(RED, 24);
//analogWrite(GREEN, 3);
//analogWrite(BLUE, 1);
//delay(1000);
  if(Serial.available() > 0){
    int out=Serial.read();
    analogWrite(RED, int(2.4*out));
    analogWrite(GREEN, int(0.3*out));
    analogWrite(BLUE, int(0.1*out));
  }
delay(random(100));
}
