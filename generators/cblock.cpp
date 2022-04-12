#include "cblock.h"

Cblock::Cblock(int order)
{
  p_order = order;
  p_xmin  = p_ymin = -1000.0;
  p_xmax  = p_ymax =  1000.0;
}

Cblock::Cblock(void)
{
  Cblock(1);
}

Cblock::~Cblock(void)
{
}

/* Set the parameters for the transfer function

  Y(s)    b0s + b1
 ----- = --------------
  U(s)    a0s + a1

  a = [a0,a1], b = [b0,b1]
 */
int Cblock::setparams(double *a,double *b)
{
  p_A[0] = -a[1]/a[0];
  p_B[0] = (b[1] - a[1]*b[0])/a[0];
  p_C[0] = 1.0;
  p_D[0] = b[0]/a[0];
}

int Cblock::setxlimits(double xmin,double xmax)
{
  p_xmin = xmin;
  p_xmax = xmax;

  return 0;
}

int Cblock::setylimits(double ymin,double ymax)
{
  p_ymin = ymin;
  p_ymax = ymax;

  return 0;
}

double Cblock::getderivative(double x, double u)
{
  return p_A[0]*x + p_B[0]*u;
}

double Cblock::updatestate(double u,double dt,double xmin,double xmax,DeltaModeStage stage)
{
  double x,dx_dt;

  if(stage == PREDICTOR) {
    p_dxdt[0] = getderivative(p_x[0],u);
    x = p_x[0] + dt*p_dxdt[0];
    x = std::max(xmin,std::min(x,xmax));
    p_xhat[0] = x;
  }

  if(stage == CORRECTOR) {
    dx_dt = getderivative(p_xhat[0],u);
    x     = p_x[0] + 0.5*dt*(p_dxdt[0] + dx_dt);
    x = std::max(xmin,std::min(x,xmax));
    p_x[0] = x;
  }

  return x;
}

double Cblock::updatestate(double u,double dt,DeltaModeStage stage)
{
  double x;

  x = updatestate(u,dt,p_xmin,p_xmax,stage);

  return x;
}

double Cblock::getoutput(double u,double dt,double xmin,double xmax,double ymin,double ymax,DeltaModeStage stage)
{
  double x,y;

  x = updatestate(u,dt,xmin,xmax,stage);
  
  y = p_C[0]*x + p_D[0]*u;

  y = std::max(ymin,std::min(y,ymax));

  return y;
}

double Cblock::getoutput(double u,double dt,DeltaModeStage stage)
{
  double y;

  y = getoutput(u,dt,p_xmin,p_xmax,p_ymin,p_ymax,stage);

  return y;
}

void Cblock::init(double u, double y)
{
  double uout;
  if(fabs(p_A[0]) < 1e-10) p_x[0] = y;
  else {
    uout = y/(p_D[0] - p_C[0]*p_B[0]/p_A[0]);
    p_x[0] = -p_B[0]/p_A[0]*uout;
  }
}

const double Cblock::getstate(DeltaModeStage stage)
{
  double x;
  if(stage == PREDICTOR) x = p_xhat[0];
  if(stage == CORRECTOR) x = p_x[0];
  return x;
}

// ------------------------------------
// PI Controller
// ------------------------------------

PIControl::PIControl(void)
{
  PIControl(1.0,1.0);
}

PIControl::PIControl(double Kp, double Ki)
{
  setconstants(Kp,Ki);
}

PIControl::PIControl(double Kp, double Ki, double xmin, double xmax, double ymin, double ymax)
{
  PIControl(Kp,Ki);

  setxlimits(xmin,xmax);
  setylimits(ymin,ymax);
}

int PIControl::setconstants(double Kp, double Ki)
{
  double a[2],b[2];

  // Parameters for state-space representation
  // Transfer funtion is
  // Kp + Ki/s => (sKp + Ki)/s
  // In standard form, this will be
  // (sKp + Ki)/(s + 0)

  b[0] = Kp;
  b[1] = Ki;
  a[0] = 1.0;
  a[1] = 0.0;

  setparams(a,b);
}
int PIControl::setconstants(double Kp, double Ki,double xmin,double xmax,double ymin,double ymax)
{
  setconstants(Kp,Ki);
  setxlimits(xmin,xmax);
  setylimits(ymin,ymax);
}

// ------------------------------------
// Filter
// ------------------------------------

Filter::Filter(void)
{
  Filter(1.0); // default time-constant is 1.0
}

Filter::Filter(double T)
{
  setconstants(T);
}

Filter::Filter(double T, double xmin, double xmax, double ymin, double ymax)
{
  Filter(Ts);

  setxlimits(xmin,xmax);
  setylimits(ymin,ymax);
}

int Filter::setconstants(double T)
{
  double a[2],b[2];

  // Parameters for state-space representation
  // Transfer funtion is
  // 1/(1 + sT)
  // In standard form, this will be
  // (s*0 + 1)/(sT + 1)

  b[0] = 0;
  b[1] = 1;
  a[0] = T;
  a[1] = 1;

  setparams(a,b);
}

int Filter::setconstants(double T,double xmin,double xmax,double ymin,double ymax)
{
  setconstants(T);
  setxlimits(xmin,xmax);
  setylimits(ymin,ymax);
}
