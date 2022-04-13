/**
 * @file cblock.h
 * @brief Header file defining the public API for linear control block
 *
 */

#ifndef CBLOCK_H
#define CBLOCK_H

#include <math.h>
#include <algorithm>

/*
  DeltaModeStage - Stage in Delta mode
*/
typedef enum {
  PREDICTOR, // Predictor update
  CORRECTOR  // Corrector update
}DeltaModeStage;


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

  State-update:
    PREDICTOR (Forward Euler):

    \hat{x}_{n+1} = x_{n} + dt*dx_dt(x_{n},u_{n})

    \hat{x}_{n+1} = max(xmin,min(\hat{x}_{n+1},xmax)

    y_{n+1} = C\hat{x}_{n+1} + Du_{n}

    y_{n+1} = max(ymin,min(y_{n+1},ymax)

    CORRECTOR (Trapezoidal):

    x_{n+1} = x_{n} + 0.5*dt*(dx_dt(x_{n}) + dx_dt(\hat{x}_{n+1},u_{n+1})) 

    x_{n+1} = max(xmin,min(x_{n+1},xmax)

    y_{n+1} = Cx_{n+1} + Du_{n+1}

    y_{n+1} = max(ymin,min(y_{n+1},ymax)

    Note here that GridLab-D does a network solve after every predictor/corrector
    call. So, during the corrector stage the input u is updated (u_{n+1})
*/ 
class Cblock
{
 protected:
  double p_A[1]; /* A */
  double p_B[1]; /* B */
  double p_C[1]; /* C */
  double p_D[1]; /* D */

  double p_dxdt[1];   /* State derivative */
  double p_xhat[1];   /* Predictor stage x */

  // p_order is kept for future extensions if and
  // when the order of the transfer function > 1
  int    p_order; /* order of the control block */

  double p_xmax,p_xmin; /* Max./Min. limits on state X */
  double p_ymax,p_ymin; /* Max./Min. limits on output Y */

  // Method for updating state x
  double updatestate(double u, double dt,double xmin, double xmax, DeltaModeStage stage);

  // Method for updating state x
  double updatestate(double u, double dt,DeltaModeStage stage);

  // Method for calculating the derivative xdot
  double getderivative(double x,double u);

 public:
  Cblock();

  // This is made public so that it can be accessed via
  // PADDR() method
  double x[1];      /* State variable x */

  // Method for setting state-space model coefficients
  void setcoeffs(double *a,double *b);

  // Method for setting x limits
  void setxlimits(double xmin, double xmax);

  // Method for setting y limits
  void setylimits(double ymin, double ymax);

  // Method for initializing state x given input u and output y
  void init(double u, double y);

  // Method for getting the output
  double getoutput(double u,double dt,DeltaModeStage stage);

  // Method for getting the ouput with limits specified
  double getoutput(double u,double dt,double xmin, double xmax, double ymin, double ymax, DeltaModeStage stage);

  // Method for getting state x
  double getstate(DeltaModeStage stage);

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

  // Methods for setting parameters
  void setparams(double Kp, double Ki);
  void setparams(double Kp, double Ki,double xmin,double xmax,double ymin,double ymax);
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
  Filter(double T, double xmin, double xmax,double ymin, double ymax);

  // Methods for setting parameters
  void setparams(double T);
  void setparams(double T,double xmin,double xmax,double ymin,double ymax);
};

#endif
