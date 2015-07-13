export void generate_triangles_and_velocities (uniform REAL lo[3], uniform REAL hi[3], uniform int nt, uniform REAL * uniform t[3][3], uniform REAL * uniform v[3], uniform unsigned int tid[], uniform unsigned int pid[])
{
  struct RNGState state;

  int seed = 0;

  while (rdrand (&seed) == false);

  seed_rng (&state, seed);

  uniform REAL scale = 1.0/pow(nt,1./3.);

  foreach (i = 0 ... nt)
  {
    REAL xrnd = frandom(&state);
    REAL yrnd = frandom(&state);
    REAL zrnd = frandom(&state);

    REAL xmul = frandom(&state)*scale;
    REAL ymul = frandom(&state)*scale;
    REAL zmul = frandom(&state)*scale;

    t[0][0][i] = xrnd + xmul*(hi[0]-lo[0]);
    t[0][1][i] = yrnd + ymul*(hi[1]-lo[1]);
    t[0][2][i] = zrnd + zmul*(hi[2]-lo[2]);

    xrnd = frandom(&state);
    yrnd = frandom(&state);
    zrnd = frandom(&state);

    xmul = frandom(&state)*scale;
    ymul = frandom(&state)*scale;
    zmul = frandom(&state)*scale;

    t[1][0][i] = xrnd + xmul*(hi[0]-lo[0]);
    t[1][1][i] = yrnd + ymul*(hi[1]-lo[1]);
    t[1][2][i] = zrnd + zmul*(hi[2]-lo[2]);

    xrnd = frandom(&state);
    yrnd = frandom(&state);
    zrnd = frandom(&state);

    xmul = frandom(&state)*scale;
    ymul = frandom(&state)*scale;
    zmul = frandom(&state)*scale;

    t[2][0][i] = xrnd + xmul*(hi[0]-lo[0]);
    t[2][1][i] = yrnd + ymul*(hi[1]-lo[1]);
    t[2][2][i] = zrnd + zmul*(hi[2]-lo[2]);

    v[0][i] = 0.5-frandom(&state);
    v[1][i] = 0.5-frandom(&state);
    v[2][i] = 0.5-frandom(&state);

    tid[i] = i;
    pid[i] = i;
  }
}

export void integrate_triangles (uniform REAL step, uniform REAL lo[3], uniform REAL hi[3],
     uniform int nt, uniform REAL * uniform t[3][3], uniform REAL * uniform v[3])
{
  foreach (i = 0 ... nt)
  {
    t[0][0][i] += step*v[0][i];
    t[0][1][i] += step*v[1][i];
    t[0][2][i] += step*v[2][i];
    t[1][0][i] += step*v[0][i];
    t[1][1][i] += step*v[1][i];
    t[1][2][i] += step*v[2][i];
    t[2][0][i] += step*v[0][i];
    t[2][1][i] += step*v[1][i];
    t[2][2][i] += step*v[2][i];

    if (t[0][0][i] < lo[0]) v[0][i] *= -1;
    if (t[0][1][i] < lo[1]) v[1][i] *= -1;
    if (t[0][2][i] < lo[2]) v[2][i] *= -1;
    if (t[0][0][i] > hi[0]) v[0][i] *= -1;
    if (t[0][1][i] > hi[1]) v[1][i] *= -1;
    if (t[0][2][i] > hi[2]) v[2][i] *= -1;

    if (t[1][0][i] < lo[0]) v[0][i] *= -1;
    if (t[1][1][i] < lo[1]) v[1][i] *= -1;
    if (t[1][2][i] < lo[2]) v[2][i] *= -1;
    if (t[1][0][i] > hi[0]) v[0][i] *= -1;
    if (t[1][1][i] > hi[1]) v[1][i] *= -1;
    if (t[1][2][i] > hi[2]) v[2][i] *= -1;

    if (t[2][0][i] < lo[0]) v[0][i] *= -1;
    if (t[2][1][i] < lo[1]) v[1][i] *= -1;
    if (t[2][2][i] < lo[2]) v[2][i] *= -1;
    if (t[2][0][i] > hi[0]) v[0][i] *= -1;
    if (t[2][1][i] > hi[1]) v[1][i] *= -1;
    if (t[2][2][i] > hi[2]) v[2][i] *= -1;
  }

}
