#ifndef SIMULATOR
#include<Wire.h>
#include<LiquidCrystal_I2C.h>
#endif

LiquidCrystal_I2C lcd(0x27, 16, 2);

#ifndef ID
#define ID 0
#endif

#define VERSION 7

#define POSITION_BUFSIZE 5

struct LinearActuator {
  int R_PWM;
  int R_EN;
  int L_PWM;
  int L_EN;
  int POS_MIN;
  int POS_MAX;
  int POS_PIN;
  double default_speed;
  double position_buffer[POSITION_BUFSIZE];
  int position_buffer_i;
  double position;
  double target;
  bool enable;
  unsigned long lastTime;
  double cError;
  double lastError;
  int last_order;
  unsigned long last_order_time;
  double pos_start;
  double speed_correction;
  bool debug;

};

#define MAX_PWM 245

double Kp=2000;
double Ki=50;
double Kd=0;

#define SPEED_CONTROL 0
//#define PLOT_PID 

#define MAX_RAMPING 800 // correspond to max variation on order per sec.

//#define DEBUG_RAMPING

LinearActuator actuators[3]={ {3,2,5,4,370,960,A0,80} ,
                              {6,7,9,8,370,960,A1,30} ,
                              {10,12,11,13,370,960,A2,80} };

#ifdef  SIMULATOR
#define debug(...) if (legId==0) {fprintf(stderr,"[%d]",legId);fprintf(stderr,__VA_ARGS__);}
#else
#define debug(...)
#endif

void setup(LinearActuator *a){
   pinMode(a->R_PWM,OUTPUT); 
   pinMode(a->L_PWM,OUTPUT); 
   pinMode(a->R_EN ,OUTPUT); 
   pinMode(a->L_EN ,OUTPUT); 
   digitalWrite(a->R_PWM,LOW);
   digitalWrite(a->L_PWM,LOW);
   digitalWrite(a->R_EN ,LOW);
   digitalWrite(a->L_EN ,LOW);
   a->lastTime = millis();
   a->enable = false;
   a->cError=0;
   a->lastError=0;
   a->last_order=0;
   a->speed_correction=0.0;
   a->debug=false;
   a->last_order_time=millis();
   for(int i=0;i<POSITION_BUFSIZE;++i)
     a->position_buffer[i]=0;
   a->position_buffer_i=0;
}

int legId;

void setup() {
  int i;
  // put your setup code here, to run once:
  Serial.begin(115200); 
  for(i=0;i<3;++i){
    setup(actuators+i);
  }  
  lcd.begin();
  lcd.backlight();
  lcd.print("Hello, world!");

  legId=ID;
  // read the id using pin A3
  //int v=analogRead(A3);
  //if (v<100) legId=0;
  //else if (v>1000) legId=3;
}

int i=0;

int qsort_dbl(const void *a,const void *b){
  if ((*((double *)a)) < (*((double *)b)))
    return -1;
  if ((*((double *)a)) < (*((double *)b)))
    return 1;
  return 0;
}

void updatePosition(LinearActuator *p){
  //if (p->enable==false) return ;
  int value=analogRead(p->POS_PIN);
  if (value<p->POS_MIN) value=p->POS_MIN;
  if (value>p->POS_MAX) value=p->POS_MAX;
  p->position_buffer[p->position_buffer_i]=value;
  p->position_buffer_i=(p->position_buffer_i+1) % POSITION_BUFSIZE;
  
  double t[POSITION_BUFSIZE]; // kmean
  for(int j=0;j<POSITION_BUFSIZE;++j)
    t[j]=p->position_buffer[j];
  qsort(t,POSITION_BUFSIZE,sizeof(t[0]),qsort_dbl);
  
  p->position = 200 *(1.0- ((double) t[POSITION_BUFSIZE/2] - (double)p->POS_MIN) / (double) (p->POS_MAX-p->POS_MIN));
}

bool computeOrder(LinearActuator *p){
  updatePosition(p);
  if (p->debug) return false;
  if (p->enable==false){
    p->speed_correction=0.0;
    p->last_order=0;
    digitalWrite(p->R_EN ,LOW);
    digitalWrite(p->L_EN ,LOW);
    return false;
  }
  unsigned long current=millis();
  double elapsed = (double)(current-p->lastTime)/1000.0;


      
  double error = (p->target - p->position)/200.0; // error is between -1 and 1

  if (fabs(error)<(2.0/200.0)){ // we arrive at target with at most 2mm of error, stop actuator
    p->enable=false;
    p->speed_correction=0.0;
    p->last_order=0;
    digitalWrite(p->R_EN ,LOW);
    digitalWrite(p->L_EN ,LOW);
    return false;
  }

  if (elapsed<(1.0/20.0)) return false;  // pid run at 40Hz
  
  double out;
  
  if (error>0) out= MAX_PWM;
  else out=-MAX_PWM;
  
  /*
  p->cError += error ;
  double rError = (error - p->lastError);

  double out=Kp * error  + Ki * p->cError + Kd * rError;
  Serial.print("D error:");Serial.print(error);
  Serial.print(" cError:");Serial.print(p->cError);
  Serial.print(" rError:");Serial.print(rError);
  Serial.print(" => ");Serial.println(out);
  p->lastError=error;

  
  //digitalWrite(p->R_EN ,HIGH);
  //digitalWrite(p->L_EN ,HIGH);

  */
  p->lastTime = current;
    
  double dts=(double )(current - p->last_order_time) / 1000.0;
  
  if (out>0){
    if (out>MAX_PWM) out=MAX_PWM;
    if (out<10) out=0;
    if (p->last_order<0){
      out=100;
    }
    else {
      double delta = (double )(out - p->last_order) / dts;
      if (delta>MAX_RAMPING)
	out = p->last_order + MAX_RAMPING * dts;
    }
    //analogWrite(p->R_PWM,out);
    //analogWrite(p->L_PWM,0);    
    
  } else {
    if (out<-MAX_PWM) out=-MAX_PWM;
    if (out>-10) out=0;
    if (p->last_order>0) {
      out=-100;
    } else {
      double delta = (double )(out - p->last_order) / dts;
      if (delta<-MAX_RAMPING)
	out=p->last_order - MAX_RAMPING * dts;
    }
    //analogWrite(p->R_PWM,0);
    //analogWrite(p->L_PWM,-out);
  }
  p->last_order=out;
  p->last_order_time=current;
  return true;
}

void do_debug(void ){
  for(int i=0;i<3;++i){
    if (actuators[i].debug){
      digitalWrite(actuators[i].R_EN ,HIGH);
      digitalWrite(actuators[i].L_EN ,HIGH);
      if (actuators[i].target>500){
	analogWrite(actuators[i].R_PWM,(int)(actuators[i].target-500));
	analogWrite(actuators[i].L_PWM,0);	
      } else {
	analogWrite(actuators[i].R_PWM,0);
	analogWrite(actuators[i].L_PWM,(int)(actuators[i].target));
      }      
    }
  }
}



void synchronizeOrders(){
  int i,j=0;
  // compute error:
  double prog[3]={-1,-1,-1};
  double tmp;
  int min_s=0;
  bool initial=false;
  debug("SYNCHRONIZE: \n");
  for(i=0;i<3;++i)
    if (actuators[i].enable){
      if (actuators[i].speed_correction==0){
	// actuator doesnot move yet
	prog[i]=fabs(actuators[i].target-actuators[i].pos_start)/actuators[i].default_speed;
	initial=true;
      }
      else
	// adjust speed according to actuators progression
	if (fabs(actuators[i].target-actuators[i].pos_start)<0.01) prog[i]=-1;
	else{
	  prog[i]=fabs(actuators[i].position - actuators[i].pos_start)/fabs(actuators[i].target-actuators[i].pos_start);
	  j+=1;
	}
      //debug("progression of %d is %lf (%d)\n",i,prog[i],j);
    }
  
  for(i=0;i<3;++i){
    debug("prog[%d]=%lf (%lf => %lf => %lf)\n",i,prog[i],actuators[i].pos_start,actuators[i].position,actuators[i].target);
  }
  
  if (j>1){
    // compute speed correction:
    if (initial){ // establish speed correction based on theoretical actuator speed:
      //prog[i] contains estimated time to make order at full speed
      debug("initial speed computation:\n");
      int slowest=0;
      for(i=1;i<3;++i){
	if ((prog[slowest]<0) || ( (prog[i]>prog[slowest]) && (prog[i]>0)))
	  slowest=i;
      }
      debug("slowest is %d \n",slowest);
      // set full speed at 1.0 for the slowest actuator:
      actuators[slowest].speed_correction=1.0;
      // use a ratio for others:
      for(i=0;i<3;++i)
	if ((i!=slowest) && (prog[i]>0))
	  actuators[i].speed_correction = prog[i]/prog[slowest];
      for(i=0;i<3;++i) {debug("%d: %lf => %lf\n",i,prog[i],actuators[i].speed_correction);}
    }
    else{
      debug("progression based speed correction:\n");
      min_s=0;
      float mm=0;
      for(i=0;i<3;++i){
	if ((prog[i]>=0) && ((prog[i]<prog[min_s])||(prog[min_s]<0)))
	  min_s=i;
	if (prog[i]>mm) mm=prog[i];
      }
      debug("slowest is %d with %lf\n",min_s,prog[min_s]);
    
      if (mm>0.01){ // Do not change inital speed correction until at leat 10% have been made so we have data to make prediction 
	// Set speed correction to 1 for the slowest
	actuators[min_s].speed_correction=1.0;
	// Asserv other to the slowest speed;
	for(i=0;i<3;++i){
	  if ((i!=min_s) && (prog[i]>=0)){
	    tmp=1.0/(prog[i] * (1.0-prog[min_s])/(prog[min_s] * (1.0 - prog[i])));
	    actuators[i].speed_correction=SPEED_CONTROL*actuators[i].speed_correction + (1.0-SPEED_CONTROL)*tmp;
	  }
	}
      }
    }
    
    // update last_order using speed_correction:
    for(i=0;i<3;++i)
      if (actuators[i].enable){
	debug("SC: %i: %lf\n",i,actuators[i].speed_correction);
	if (actuators[i].speed_correction<0) actuators[i].speed_correction=0.1;
	actuators[i].last_order*=actuators[i].speed_correction;
	if (actuators[i].last_order>MAX_PWM) actuators[i].last_order=MAX_PWM;
	if (actuators[i].last_order<-MAX_PWM) actuators[i].last_order=-MAX_PWM;
	//debug("=> %d\n",actuators[i].last_order);
      }    
  } else { // reset speed correction as there is only one active actuator and let last_order as computed
    for(i=0;i<3;++i)
      actuators[i].speed_correction=1.0;    
  }
  
  // send order:
  for(i=0;i<3;++i){    
    if (actuators[i].enable){      
      digitalWrite(actuators[i].R_EN ,HIGH);
      digitalWrite(actuators[i].L_EN ,HIGH);
      if (actuators[i].last_order>0){
	analogWrite(actuators[i].R_PWM,actuators[i].last_order);
	analogWrite(actuators[i].L_PWM,0);
      } else {
	analogWrite(actuators[i].R_PWM,0);
	analogWrite(actuators[i].L_PWM,-actuators[i].last_order);
      }
    }
  }
}


char buffer[100];
int r=0;
int v;

void update_target(char *order){
  LinearActuator *p;
  if (order[0]=='A')      { p=&(actuators[0]); p->debug=false;}
  else if (order[0]=='a') { p=&(actuators[0]); p->debug=true;}
  else if (order[0]=='B') { p=&(actuators[1]); p->debug=false;}
  else if (order[0]=='b') { p=&(actuators[1]); p->debug=true;}
  else if (order[0]=='C') { p=&(actuators[2]); p->debug=false;}
  else if (order[0]=='b') { p=&(actuators[2]); p->debug=true;}
  else return;
  p->target=0;
  if (order[1]=='_') {
    //Serial.println("disabled");
    p->debug=false;
    p->enable=false;
    p->speed_correction=0.0;
    p->last_order=0;
    return;
  }
  
  int i=1;
  double deci=1;
  bool d=false;
  while((order[i]!='\n') && (order[i]!='\0') && (order[i]!='#')){    
    if (order[i]=='.'){
      d=true;
    } else {
      p->target=p->target*10;
      p->target+=order[i]-'0';
      if (d) deci*=10;
    }
    i+=1;
  }
  p->target = p->target / deci;
  Serial.print("D new target is: (");Serial.print(millis());Serial.print("): ");Serial.println(p->target);
  p->speed_correction=0.0;
  p->pos_start = p->position;
  if (p->enable==false){
    p->enable=true;
    p->lastTime=millis();
    p->cError=0;
    p->lastError=0;
    p->last_order=0;
    p->last_order_time=millis();
    p->pos_start = p->position;
  }
  if (order[i]!='\0')
    update_target(order+i+1);
}

long lcdfreq=0;

int lcdi=0;

long statusfreq=0;
int statusBlink=0;

void loop() {
  if ((millis()-lcdfreq)>1000){
    lcd.clear();
    lcd.print("leg number: ");lcd.print(legId);lcd.print("(");lcd.print(VERSION);lcd.print(")");
    lcd.setCursor(0,1);
    if (actuators[lcdi].enable)
      {lcd.print((char)('A'+lcdi));lcd.print(":");lcd.print(actuators[lcdi].position); lcd.print("/"); lcd.print(actuators[lcdi].target);}
      else{lcd.print((char)('a'+lcdi)); lcd.print(": disabled");}
    lcdfreq=millis();
    lcdi=(lcdi+1)%3;
  }

  if ((millis()-statusfreq)>200){ // send status every 200ms
    /*
    if (statusBlink==0)
      lcd.noBacklight();
    else if (statusBlink==10)
      lcd.backlight();
    statusBlink=(statusBlink+1)%20;
    */
      
    Serial.print("S ");
    for(int i=0;i<3;++i){
      Serial.print(actuators[i].enable?1:0);
      Serial.print(';');
      Serial.print(actuators[i].position);
      Serial.print(';');
      Serial.print(actuators[i].target);
      Serial.print(';');
      Serial.print(actuators[i].last_order);
      Serial.print(' ');
    }
    Serial.println();
    statusfreq=millis();
    /*    Serial.print("D speed:");
    for(int i=0;i<3;++i){
      Serial.print(actuators[i].speed_correction);
      Serial.print(';');
    }
    Serial.println();
    */
  }
  v=Serial.read();
#ifdef SIMULATOR
  //if (v!=0) fprintf(stderr,"verrin read %c \n",v);
#endif
  if (v=='?'){
    Serial.print("id[");
    Serial.print(legId);
    Serial.println("]");
  }
  if ((r==0 && ((v>='A' && v<='C') || ((v>='a')&&(v<='c')) )) || (r>0 && (isalpha(v) || isdigit(v) || v=='.' || v=='_' || v=='#'))){
    //    if (((v>='A') && (v<='C')) || ((v>='a')&&(v<='c')))
    //      r=0;
    buffer[r]=v;
    r+=1;
  } 
  if (v=='\n'){
    buffer[r]='\0';
    update_target(buffer);
    r=0;
  }

  bool todo=false;
  todo= computeOrder(actuators+0) || todo;
  todo= computeOrder(actuators+1) || todo;
  todo= computeOrder(actuators+2) || todo;
  if (todo)
    synchronizeOrders();
  do_debug();
}
