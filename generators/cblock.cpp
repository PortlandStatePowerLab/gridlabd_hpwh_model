#include "cblock.h"

Cblock::Cblock(int order)
{
  p_order = order;
  p_hasxlimiter = false;
  p_xatmax = p_xatmin = false;
  p_hasylimiter = false;
  p_updatestatecalled = false;
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

  p_hasxlimiter = true;

  return 0;
}

int Cblock::setylimits(double ymin,double ymax)
{
  p_ymin = ymin;
  p_ymax = ymax;

  p_hasylimiter = true;

  return 0;
}

double Cblock::getderivative(double x, double u)
{
  return p_A[0]*x + p_B[0]*u;
}

double Cblock::getoutput(double x, double u)
{
  double y;
  y = p_C[0]*x + p_D[0]*u;
  if(p_hasylimiter) {
    if(y < p_ymin) y = p_ymin;
    else if(y > p_ymax) y = p_ymax;
  }
  return y;
}

double Cblock::updatestate(double u,double dt)
{
  double x_n = p_x[0]; // State at time-step n
  double x_n1;      // State at time-step n+1
  double px_n1;     // predicted state with Euler approximation
  double dpx_dtn1;  // Predicted derivate for state px_n1
  double dx_dtn;     // Derivative at time-step n

  dx_dtn = getderivative(x_n,u);

  // Approximation for x_n1
  px_n1 = x_n + dt*dx_dtn;

  // Derivative
  dpx_dtn1 = getderivative(px_n1,u);
  dpx_dtn1 += dx_dtn;

  if(p_xatmin) {
    if(dpx_dtn1 > 0) { // Release
      x_n1 = x_n + 0.5*dt*dpx_dtn1;
    } else {
      x_n1 = p_xmin;
    }
  } else if(p_xatmax) {
    if(dpx_dtn1 < 0) { // Release
      x_n1 = x_n + 0.5*dt*dpx_dtn1;
    } else {
      x_n1 = p_xmax;
    }
  } else {
    x_n1 = x_n + 0.5*dt*dpx_dtn1;
    if(x_n1 > p_xmax) {
      x_n1 = p_xmax;
      p_xatmax = true;
    } else if(x_n1 < p_xmin) {
      x_n1 = p_xmin;
      p_xatmin = p_xmin;
    }
  }

  p_x[0] = x_n1; // Update state
    
  return x_n1;
}

int Cblock::init(double u, double y)
{
  double uout = y/(p_D[0] - p_C[0]*p_B[0]/p_A[0]);
  p_x[0] = -p_B[0]/p_A[0]*uout;
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
  // 1/(1 + sT) => (sKp + Ki)/s
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
