// Simple trial wave function and accompaying kinetic energy.
// Potential energy is without coloumb term.
//----------------------------------------------------------------------------

inline double acceptAmp1(double *params, double omega, mat &r, vec &delta,
                         int num) {
  double alpha = params[0];
  double expo = dot(delta, delta + 2 * r.col(num));
  return exp(-alpha * omega * expo);
}

inline double localKinetic1(double *params, double omega, mat &r) {
  double R2 = accu(mat(r % r)); // r%r means elementwise multiplication
  double alpha = params[0];
  return -0.5 * omega * omega * R2 * alpha * alpha + 3 * alpha * omega;
}

inline double localPotential1(double *params, double omega, mat &r) {
  double R2 = accu(mat(r % r)); // r%r means elementwise multiplication
  return 0.5 * omega * omega * R2;
}
//----------------------------------------------------------------------------

// Trial wave function with Jastrow factor and accompaying kinetic energy.
// Potential energy is with coloumb term.
//----------------------------------------------------------------------------
inline double acceptAmp2(double *params, double omega, mat &r, vec &delta,
                         int num) {
  int other;
  if (num == 0)
    other = 1;
  else
    other = 0;

  double alpha = params[0];
  double beta = params[1];

  double R1 = -alpha * omega * dot(delta, delta + 2 * r.col(num));
  double R2 = norm(r.col(num) + delta - r.col(other));
  double R3 = norm(r.col(num) - r.col(other));

  return exp(R1 + R2 / (1 + beta * R2) - R3 / (1 + beta * R3));
}

inline double localKinetic2(double *params, double omega, mat &r) {
  double alpha = params[0];
  double beta = params[1];

  double R12 = norm(r.col(0) - r.col(1));
  double BR = 1 + beta * R12;
  return localKinetic1(params, omega, r) +
         1 / (2 * BR * BR) *
             (alpha * omega * R12 - 1 / (2 * BR * BR) - 2 / R12 +
              2 * beta / BR);
}

inline double localPotential2(double *params, double omega, mat &r) {
  double R2 = accu(mat(r % r)); // r%r means elementwise multiplication
  return 0.5 * omega * omega * R2 + 1 / norm(r.col(0) - r.col(1));
}
//----------------------------------------------------------------------------
