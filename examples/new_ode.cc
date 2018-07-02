typedef std::array<double, dimX> X_type;
typedef std::array<double, 2*dimX> Y_type;

auto Next = [](Y_type &y, U_type &u, double tau, OdeSolver solver) -> void {
  auto rhs =[](Y_type &yy, const Y_type &y, U_type &u) -> void {
    /* find the distrubance for the given state */
    X_type x;
    for (int i=0; i<dimX; i++)
      x[i] = y[i];
    X_type w = disturbance(x);

    /* coupled system + growth bound ode */
    yy[0] = u[0]*std::cos(y[2]);
    yy[1] = u[0]*std::sin(y[2]);
    yy[2] = u[1];
    yy[3] = y[5]*std::abs(u[0]) + w[0];
    yy[4] = y[5]*std::abs(u[0]) + w[1];
    yy[5] = w[2];
  };
  solver(rhs, y, u);
};
