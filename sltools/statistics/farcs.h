#ifndef FARCS_H
#define FARCS_H


double f_arcs_I_arg(double zs, double params[], CrossSection A){
  double zl = params[0];
  double M = params[1];
  double conc = params[2];
  double m_lim = params[3];
  double cosmological_params[] = {params[4],params[5],params[6]};
  //printf("zs=%E,zl=%E,M=%E,Conc=%E,m_lim=%E\n,Om=%E,Ok=%E,OL=%E\n",zs,zl,M,conc,m_lim,cosmological_params[0],cosmological_params[1],cosmological_params[2] );
  //printf("%E\n",DimCrossSection(zl,zs,M,conc,cosmological_params,A));
  //printf("%E\n",n(zs,m_lim));

  return DimCrossSection(zl,zs,M,conc,cosmological_params,A)*n(zs,m_lim);
}

double f_arcs_I(double zl, double M, double conc, double m_lim, double cosmological_params[], CrossSection A){
  double params[] = {zl, M, conc, m_lim, cosmological_params[0], cosmological_params[1], cosmological_params[2]};

  double INT = 0.0;
  double D = 0.0;

  double z_bean = 0.05;
  double zf = zl+z_bean;

  for (int i=1;i<=400;i++){
    D = IntGaulegSub10(f_arcs_I_arg, params, A, zl, zf);
    INT+= D;
    if (0.0001*INT/double(i)>D) break;
    zl += z_bean;
    zf += z_bean;
  }
  return INT;
  //return IntGaulegDef(f_arcs_I_arg , params, A, zl, 1000, 100);
}


double f_arcs_II_arg(double zs, double params[], CrossSection A){
  double zl = params[0];
  double M = params[1];
  double conc = params[2];
  double cosmological_params[] = {params[4],params[5],params[6]};

  double ks = kappa_s(zl, zs, M, conc, cosmological_params);
  double m_lim = 1.25*log10(A.mag_med(ks)) + params[3];
  //printf("zs=%E,zl=%E,M=%E,Conc=%E,m_lim=%E\n,Om=%E,Ok=%E,OL=%E\n",zs,zl,M,conc,m_lim,cosmological_params[0],cosmological_params[1],cosmological_params[2] );
  //printf("%E\n",DimCrossSection(zl,zs,M,conc,cosmological_params,A));
  //printf("%E\n",n(zs,m_lim));
  return DimCrossSection(zl,zs,M,conc,cosmological_params,A)*n(zs,m_lim);
}


double f_arcs_II(double zl, double M, double conc, double m_lim, double cosmological_params[], CrossSection A){
  double params[] = {zl, M, conc, m_lim, cosmological_params[0], cosmological_params[1], cosmological_params[2]};

  double INT = 0.0;
  double D = 0.0;

  double z_bean = 0.05;
  double zf = zl+z_bean;

  for (int i=1;i<=400;i++){
    D = IntGaulegSub10(f_arcs_II_arg, params, A, zl, zf);
    INT+= D;
    if (0.0001*INT/double(i)>D) break;
    zl += z_bean;
    zf += z_bean;
  }
  return INT;
  //return IntGaulegDef(f_arcs_II_arg , params, A, zl, 20, 1000);
}

double f_arcs_III_arg_mu(double mu, double params[], CrossSection A){
  double zl = params[0];
  double zs = params[1];
  double M = params[2];
  double conc = params[3];
  double cosmological_params[] = {params[5],params[6],params[7]};

  double ks = kappa_s(zl, zs, M, conc, cosmological_params);
  double mu_med = A.mag_med(ks);
  double m_lim = 1.25*log10(mu) + params[4];
  //printf("zs=%E,zl=%E,M=%E,Conc=%E,m_lim=%E\n,Om=%E,Ok=%E,OL=%E\n",zs,zl,M,conc,m_lim,cosmological_params[0],cosmological_params[1],cosmological_params[2] );
  //printf("%E\n",DimCrossSection(zl,zs,M,conc,cosmological_params,A));
  //printf("%E\n",n(zs,m_lim));
  return DimCrossSection(zl,zs,M,conc,cosmological_params,A)*n(zs,m_lim)*(mu_med*mu_med/(2.0*mu*mu*mu));
}

double Int_f_arcs_III_arg_mu(double zl, double zs, double M, double conc, double m_lim, double cosmological_params[], CrossSection A){
  double params[] = {zl, zs, M,conc, m_lim, cosmological_params[0], cosmological_params[1], cosmological_params[2]};

  /*double INT = 0.0;
  double D = 0.0;*/

  double ks = kappa_s(zl, zs, M, conc, cosmological_params);
  double mu_min = A.mag_med(ks)/2.0;

  /*double mu_bean = 25.0;
  double mu_final = mu_min + mu_bean;

  for (int i=1;i<=10000000;i++){
    D = IntGaulegSub10(f_arcs_III_arg_mu, params, A, mu_min, mu_final);
    INT+= D;
    if (0.01*INT/double(i)>D) break;
    mu_min += mu_bean;
    mu_final += mu_bean;
  }*/

  //return INT;
  return IntGaulegDef(f_arcs_III_arg_mu , params, A,mu_min, 1000, 150);
}


double f_arcs_III_arg(double zs, double params[], CrossSection A){
  double zl = params[0];
  double M = params[1];
  double conc = params[2];
  double cosmological_params[] = {params[4],params[5],params[6]};

  //double ks = kappa_s(zl, zs, M, conc, cosmological_params);
  double m_lim = params[3];

  return Int_f_arcs_III_arg_mu(zl, zs, M, conc, m_lim, cosmological_params, A);
}


double f_arcs_III(double zl, double M, double conc, double m_lim, double cosmological_params[], CrossSection A){
  double params[] = {zl, M, conc, m_lim, cosmological_params[0], cosmological_params[1], cosmological_params[2]};

  double INT = 0.0;
  double D = 0.0;

  double z_bean = 0.5;
  double zf = zl+z_bean;

  for (int i=1;i<=40;i++){
    D = IntGaulegSub10(f_arcs_III_arg, params, A, zl, zf);
    INT+= D;
    if (0.01*INT/double(i)>D) break;
    zl += z_bean;
    zf += z_bean;
  }
  return INT;
 // return IntGaulegSub100(f_arcs_III_arg , params, A, zl, 10000);
}
#endif





































