/**
 * @file cblock.h
 * @brief Header file defining the public API for linear control block
 *
 */

#ifndef CBLOCK_H
#define CBLOCK_H

#include <math.h>

/*
  Linear control block base class - This is the base class for linear blocks. All 
  control blocks inherit from this class. It expresses a first order linear transfer
  function in observable canonical form to obtain a first order state-space 
  representation. The state x and output y can be bounded by limits.

  Linear control block:
  input  : u 
  output : y
  state  : x
      
                 xmax
                ----
               /             ymax
              /             -----
        -------------      /
        | b0s + b1  |     /
  u ----| --------- |----------- y
        | a0s + a1  |    /
        -------------   /
             /       ---
            /        ymin
        ----
        xmin


  For a first order transfer function, 

  Y(s)    b0s + b1
  ----- = ----------
  U(s)    a0s + a1

  The equivalent state-space representation in the observable canonical form is given by

  dx_dt = Ax + Bu
    y   = Cx + Du

  here u,x,y \in R^1 

  A = -a1/a0, B = b1 - a1b0, C = 1, D = b0/a0

  Output:
   y = Cx + D*u,  ymin <= y <= ymax

  State-update: (using Modified-Euler method)

   *** x within limits ***

   if xmin < x_n < xmax
    px_{n+1} = x_{n} + dt*dx_dt_{n}
    cx_{n+1}  = x_{n} + 0.5*dt*(dx_dt(x_{n}) + dx_dt(px_{n+1})

    if x < xmin
      x_{n+1} = xmin // Snap
    elif x > xmax
      x_{n+1} = xmax // Snap
    else
      x_{n+1} = cx_{n+1}

   *** x at limit ***
   px_{n+1} = x_{n} + dt*dx_dt_{n} 
   dx_dt_{n} = dx_dt(x_{n}) + dx_dt(px_{n+1})

   if x at xmin
     if dx_dt_{n} > 0
       x_{n+1} = x_n + 0.5*dt*dx_dt_{n} // Release
     else
       x_{n+1} = xmin
   else if x at xmax
     if dx_dt_{n} < 0
       x_{n+1} = x_n + 0.5*dt*dx_dt_{n} // Release
     else
       x_{n+1} = xmax

*/ 
class Cblock
{
 protected:
  double p_A[1]; /* A */
  double p_B[1]; /* B */
  double p_C[1]; /* C */
  double p_D[1]; /* D */
  bool   p_xatmax,p_xatmin; /* Flags for max./min. limits */
  double p_x[1]; /* State variable x */
  double p_xmax,p_xmin; /* Max./Min. limits on state X */
  double p_ymax,p_ymin; /* Max./Min. limits on output Y */
  bool   p_hasxlimiter; /* Model has x limiter */
  bool   p_hasylimiter; /* Model has y limiter */
  bool   p_updatestatecalled; /* Has state been updated for this time-step */
  // p_order is kept for future extensions if and
  // when the order of the transfer function > 1
  int    p_order; /* order of the control block */

 public:
  Cblock();
  Cblock(int);
  // Method for settting state-space parameters
  int setparams(double *a,double *b);

  // Method for setting x limits
  int setxlimits(double xmin, double xmax);

  // Method for setting y limits
  int setylimits(double ymin, double ymax);

  // Method for initializing state x given input u and output y
  double init(double u, double y);

  // Method for getting the output
  double getoutput(double x,double u);

  // Method for calculating the derivative xdot
  double getderivative(double x,double u);

  // Method for updating state x
  double updatestate(double u, double dt);

  // Method for getting state x
  const double getstate();

  ~Cblock(void);
};

/*
  PI control block:
  input  : u 
  output : y
  state  : x (integrator)
      
                 xmax
                ----
               /             ymax
              /             -----
        -------------      /
        |           |     /
  u ----| Kp + Ki/s |----------- y
        |           |    /
        -------------   /
             /       ---
            /        ymin
        ----
        xmin

   Differential equation:
       dx_dt = Ki*u

   Output:
    y = x + Kp*u,  ymin <= y <= ymax

*/
class PIControl: public Cblock
{
 public:
  PIControl();
  PIControl(double Kp, double Ki);
  PIControl(double Kp, double Ki, double xmin, double xmax);
  PIControl(double Kp, double Ki, double xmin, double xmax,double ymin, double ymax);

  // Methods for setting constants
  int setconstants(double Kp, double Ki);
  int setconstants(double Kp, double Ki,double xmin,double xmax,double ymin,double ymax);
};

/*
  Filter block:
  input  : u 
  output : y
  state  : x
      
                 xmax
                ----
               /             ymax
              /             -----
        -------------      /
        |    1      |     /
  u ----| -------   |----------- y
        |  1 + sT   |    /
        -------------   /
             /       ---
            /        ymin
        ----
        xmin

   Differential equation:
       dx_dt = (u - x)/T

   Output:
    y = x,  ymin <= y <= ymax

*/
class Filter: public Cblock
{
 public:
  Filter();
  Filter(double T);
  Filter(double T, double xmin, double xmax);
  Filter(double T, double xmin, double xmax,double ymin, double ymax);

  // Methods for setting constants
  int setconstants(double T);
  int setconstants(double T,double xmin,double xmax,double ymin,double ymax);
};



#endif
