#define SIMULATOR
#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>
#include <sys/select.h>
#include <chrono>
#include <unistd.h>
#include <math.h>

using namespace std::chrono;

class LiquidCrystal_I2C{
public:
  LiquidCrystal_I2C(int,int,int){}
  void begin(){}
  void backlight(){}
  void noBacklight(){}
  void print(const char *){}
  void print(char ){}
  void println(const char *){}
  void println(char ){}
  void println(){}
  void clear(){}
  void setCursor(int,int){}
};


class SerialCom{
public: 
  void begin(int baudrate){}
  void print(const char *s){printf("%s",s);}
  void println(const char *s){printf("%s\r\n",s);}
  void print(unsigned long int i){printf("%lu",i);}
  void println(unsigned long int i){printf("%lu\r\n",i);}
  void print(int i){printf("%d",i);}
  void println(int i){printf("%d\r\n",i);}
  void print(double i){printf("%lf",i);}
  void println(double i){printf("%lf\r\n",i);}
  void print(float i){printf("%f",i);}
  void println(float i){printf("%f\r\n",i);}
  void print(char c){putchar(c);}
  void println(char c){printf("%c\r\n",c);}
  void println(){printf("\r\n");}
  char read(){
    static fd_set rfds;
    static struct timeval tv={0,1};

    FD_ZERO(&rfds);
    FD_SET(0, &rfds);
    if (select(1,&rfds,NULL,NULL,&tv)>0){
      char c;
      ::read(0,&c,1);
      return c;
    }
    return 0;    
  }
};

SerialCom Serial;

#define OUTPUT 1
#define LOW 2
#define HIGH 3
#define A0 0
#define A1 1
#define A2 2

unsigned long millis(){
  static unsigned long start=duration_cast< duration<unsigned long,std::milli> >(system_clock::now().time_since_epoch()).count();
  unsigned long ms = duration_cast< duration<unsigned long,std::milli> >(system_clock::now().time_since_epoch()).count();
  return ms-start;
}

void pinMode(int pin,int mode){}



struct motor{
  double position;
  int l_pwm;
  int r_pwm;
  bool enable;
  float current_speed;
  float max_speed;
};

struct motor motors[3]={ {0,0,0,false,0,80} , {0,0,0,false,0,40}, {0,0,0,false,0,80} };



double fabs(double a){
  if (a<0) return -a;
  return a;
}

void analogWrite(int pin,int value);
int analogRead(int pin);
void digitalWrite(int pin, int value);

#include "./verrin.ino"


void digitalWrite(int pin, int value){
  // must be hack to know if motor is enabled
  for(int i=0;i<3;++i)
    if (actuators[i].R_EN==pin) {
      if (value==LOW) motors[i].enable=false;
      else motors[i].enable=true;
    }
}

int analogRead(int pin){
  // must be hack to give potar value
  int v;
  for(int i=0;i<3;++i)
    if (actuators[i].POS_PIN==pin){
      v=(200.0-motors[i].position)/200.0*double(actuators[i].POS_MAX-actuators[i].POS_MIN)	
	+actuators[i].POS_MIN;
      //printf("%lf => %d \n",motors[i].position,v);
      return v;
    }
  return 0;
}

void analogWrite(int pin,int value){
  // must be hack to know to simulate motor action
  for(int i=0;i<3;++i){
    if (actuators[i].R_PWM==pin)
      motors[i].r_pwm = value;
    if (actuators[i].L_PWM==pin)
      motors[i].l_pwm = value;
  }
}


void move_motors(){
  static unsigned long last=millis();
  float elapse = millis()- last;
  if (elapse<5) return ;
  last=millis();
  for(int i=0;i<3;++i)
    if (motors[i].enable){
      if (motors[i].current_speed==0) motors[i].current_speed=10;
      else motors[i].current_speed=motors[i].current_speed + 2*elapse/1000.0;
      if (motors[i].current_speed>motors[i].max_speed) motors[i].current_speed=motors[i].max_speed;
      //printf("%d : %d %d %lf\n",i,motors[i].r_pwm,motors[i].l_pwm,motors[i].current_speed);
      if ((motors[i].r_pwm>30) && (motors[i].position<200))
	motors[i].position+= (motors[i].current_speed * (double)motors[i].r_pwm/(double)MAX_PWM  * (double)(elapse) ) / 1000.0;
      else if ((motors[i].l_pwm>30) && (motors[i].position>0))
	motors[i].position-=(motors[i].current_speed * (double)motors[i].l_pwm/(double)MAX_PWM * (double)(elapse) ) / 1000.0;
      // add random noise:
      motors[i].position+=(double)(rand()%50 - 25)/50.0 * 0.1; // +- 0.1mm of noise
      
    }
}



int main(int argc, char **argv){


  setvbuf(stdin,NULL,_IONBF,0);
  setvbuf(stdout,NULL,_IONBF,0);
  
  
  setup();
  legId=atoi(argv[1]);

  while(1){
    loop();
    move_motors();
  }




}



