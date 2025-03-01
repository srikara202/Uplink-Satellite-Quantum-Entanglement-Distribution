#include <iostream>
#include <cmath>
#include <functional>

//////////////////////////////////////////////
// Global/Useful Constants and Helpers
//////////////////////////////////////////////

static const double c_      = 299792458.0;            // Speed of light (m/s)
static const double ER      = 6.371e6;                // Earth's radius (m)
static const double pi      = 3.14159265358979323846; // Pi
static const double euler_e = 2.71828182845904523536; // e

double integrate(const std::function<double(double)>& f,
                 double a, double b,
                 int nSteps = 2000)
{
    // Basic check
    if(a == b) return 0.0;
    // Trapezoidal rule
    double h = (b - a)/nSteps;
    double sum = 0.5*( f(a) + f(b) );
    for(int i = 1; i < nSteps; i++){
        double x = a + i*h;
        sum += f(x);
    }
    return sum * h;
}

using std::erf;

////////////////////////////////////////////////////
// 1. Wavepacket / Photonic definitions
////////////////////////////////////////////////////

// The fundamental Gaussian wavepacket:

double psi(double x, double x0, double sigma)
{
    // (1/(2*pi*sigma^2))^(1/4)
    double prefactor = std::pow(1.0/(2.0*pi*sigma*sigma), 0.25);
    double exponent  = -0.25*std::pow((x - x0)/sigma, 2);
    return prefactor * std::exp(exponent);
}

double psi1(double x, double x0, double sigma)
{
    return psi(x, 0.0, sigma);
}

double psi2(double x, double x0, double sigma)
{
    return psi(x, x0, sigma);
}

// Probability that photon passes the gating window [c*tmin, c*tmax], i.e. integral of |psi1|^2

double Pgw1(double x0, double sigma, double tmin, double tmax)
{
    double a = c_*tmin;
    double b = c_*tmax;
    // integrand = |psi1(x)|^2
    auto integrand = [x0, sigma](double x){
        double val = psi1(x, x0, sigma);
        return val*val;
    };
    return integrate(integrand, a, b);
}

// Probability that photon passes the gating window [c*tmin, c*tmax], i.e. integral of |psi2|^2

double Pgw2(double x0, double sigma, double tmin, double tmax)
{
    double a = c_*tmin;
    double b = c_*tmax;
    auto integrand = [x0, sigma](double x){
        double val = psi2(x, x0, sigma);
        return val*val;
    };
    return integrate(integrand, a, b);
}

// Overlap integral gamma
double GammaVal(double x0, double sigma, double tmin, double tmax)
{
    double a = c_*tmin;
    double b = c_*tmax;
    auto integrand = [x0, sigma](double x){
        return psi1(x, x0, sigma)*psi2(x, x0, sigma);
    };
    double val = integrate(integrand, a, b);
    return val*val; // absolute square
}

// Fidelity due to mode mismatch: 
double Fic(double x0, double sigma, double tmin, double tmax)
{
    double pgw1 = Pgw1(x0, sigma, tmin, tmax);
    double pgw2 = Pgw2(x0, sigma, tmin, tmax);
    double gamma = GammaVal(x0, sigma, tmin, tmax);
    // Avoid zero-division:
    if(pgw1 <= 0.0 || pgw2 <= 0.0) return 0.5; 
    return 0.5 + gamma/(2.0*pgw1*pgw2);
}

////////////////////////////////////////////////////
// 2. Geometries, angles, distances, etc.
////////////////////////////////////////////////////

double zTheta(double h, double theta)
{
    return std::sqrt(h*h + 2.0*h*ER + ER*ER*std::cos(theta)*std::cos(theta))
           - ER*std::cos(theta);
}

double zAlpha(double h, double alpha)
{
    double R  = ER + h;
    double cA = std::cos(alpha);
    return std::sqrt(ER*ER + R*R - 2.0*ER*R*cA);
}

double AlphaAngle(double h, double theta)
{
    double numer = ER + zTheta(h, theta)*std::cos(theta);
    double denom = ER + h;
    double ratio = numer/denom;
    // Clip ratio to [-1,1] if needed:
    if(ratio > 1.0) ratio = 1.0;
    if(ratio < -1.0) ratio = -1.0;
    return std::acos(ratio);
}

double hTheta(double z, double theta)
{
    return std::sqrt(ER*ER + z*z + 2.0*z*ER*std::cos(theta)) - ER;
}


double alpha2(double h, double theta1, double DG)
{
    return (DG/ER) - AlphaAngle(h, theta1);
}

double theta2(double h, double theta1, double DG)
{
    double a2 = alpha2(h, theta1, DG);
    double R  = ER + h;
    double cA2 = std::cos(a2);
    double zA2 = zAlpha(h, a2);
    double ratio = (R*cA2 - ER)/zA2;
    if(ratio>1.0) ratio=1.0; 
    if(ratio<-1.0) ratio=-1.0;
    return std::acos(ratio);
}

// The "equidistant" zenith angle, theta_e, for symmetrical geometry:

double theta_e(double h, double DG)
{
    double halfDGoverER = (DG/(2.0*ER));
    double R  = ER + h;
    double cVal = std::cos(halfDGoverER);
    double top = (R*cVal - ER);
    double zA  = zAlpha(h, halfDGoverER);
    double ratio = top / zA;
    if(ratio>1.0) ratio=1.0;
    if(ratio<-1.0) ratio=-1.0;
    return std::acos(ratio);
}

////////////////////////////////////////////////////
// 3. Beam Widening and Transmission
////////////////////////////////////////////////////

// Kappa0
static const double Cw = 2.2354e-12;

double Kappa0(double theta, double lambda)
{
    double val = 1.46 * std::pow((2.0*pi/lambda),2.0) * (1.0/std::cos(theta)) * Cw;
    return std::pow(val, -3.0/5.0);
}

// wSquared
double wSquared(double w0, double z, double lambda, double K0)
{
    double term1 = w0*w0 * (1.0 + std::pow( (z*lambda)/(pi*w0*w0), 2.0));
    double term2 = 2.0 * std::pow( (lambda*z)/(pi*K0), 2.0 );
    return term1 + term2;
}

// eta_w
double eta_w(double RA, double wSq, double z, double sigma_tr)
{
    double denom = wSq + (z*1e-6)*(z*1e-6) + sigma_tr*sigma_tr;
    double arg   = -2.0*(RA*RA)/denom;
    return 1.0 - std::exp(arg);
}

////////////////////////////////////////////////////
// 4. Atmospheric Efficiency
////////////////////////////////////////////////////

static const double alpha0 = 5e-6;  // at lambda=800nm
static const double hTilde = 6600.0;


double eta_a(double h, double theta)
{
    double zt = zTheta(h, theta);
    // If zt<0, clamp it
    if(zt <= 0.0) return 1.0;
    auto integrand = [theta](double y){
        double val = hTheta(y, theta)/hTilde;
        return std::exp(-val);
    };
    double val = integrate(integrand, 0.0, zt);
    double exponent = -alpha0*val;
    return std::exp(exponent);
}

////////////////////////////////////////////////////
// 5. Overall Single-Photon Channel Efficiency
////////////////////////////////////////////////////

double eta_ph(double RA, double w0, double h, double theta, 
              double lambda, double sigma_tr, double eta_m)
{
    double K0 = Kappa0(theta, lambda);
    double z  = zTheta(h, theta);
    double wsq = wSquared(w0, z, lambda, K0);
    double eW  = eta_w(RA, wsq, z, sigma_tr);
    double eA  = eta_a(h, theta);
    return eW * eA * eta_m;
}

// For convenience, define the gating-window-limited “channel efficiencies” from A or B

double etaA(double RA, double w0, double h, double theta,
            double lambda, double sigma_tr, double eta_m,
            double x0, double sigma, double tmin, double tmax)
{
    return eta_ph(RA, w0, h, theta, lambda, sigma_tr, eta_m)
         * Pgw1(x0, sigma, tmin, tmax);
}
double etaB(double RA, double w0, double h, double theta,
            double lambda, double sigma_tr, double eta_m,
            double x0, double sigma, double tmin, double tmax)
{
    return eta_ph(RA, w0, h, theta, lambda, sigma_tr, eta_m)
         * Pgw2(x0, sigma, tmin, tmax);
}

////////////////////////////////////////////////////
// 6. Stray Photons (Day/Night)
////////////////////////////////////////////////////

// Earth albedo, solar irradiance, etc.
static const double EE = 0.3;        // Earth's Albedo
static const double IS = 4.61e27;    // Solar spectral irradiance (SI)
static const double M  = 0.14;       // Moon albedo
static const double rM = 1737.4e3;   // Moon radius (m)
static const double lME= 363300e3;   // Earth-moon dist. (m)
static const double p_ = 6.626e-34;  // Planck's const
static const double kB = 1.381e-23;  // Boltzmann's const


double ND(double RA, double thetaFOV)
{
    return EE * IS * std::pow(RA*thetaFOV, 2.0);
}

double IBB(double lambda, double T)
{
    double top = 2.0*c_/(std::pow(lambda, 4.0));
    double expo= p_*c_/(lambda*kB*T);
    double denom= std::exp(expo) - 1.0;
    return top*(1.0/denom);
}

double NN(double lambda, double T, double RA, double thetaFOV)
{
    double firstTerm = pi*IBB(lambda, T)*std::pow(RA*thetaFOV,2.0);
    double secondTerm= ND(RA,thetaFOV)* M * std::pow(rM/lME,2.0);
    return firstTerm + secondTerm;
}

double rday(double CD, double h, double theta, double eta_m,
            double RA, double thetaFOV, double DeltaLambda)
{
    double val = 0.5*eta_a(h,theta)*eta_m*ND(RA,thetaFOV)*DeltaLambda;
    return CD + val;
}

double rnight(double CD, double h, double theta, double eta_m,
              double RA, double thetaFOV, double DeltaLambda,
              double lambda, double T)
{
    double val = 0.5*eta_a(h,theta)*eta_m*NN(lambda, T, RA, thetaFOV)*DeltaLambda;
    return CD + val;
}

// Poissonian Probability of n stray photons passing in time t
//    Psp(n, r, t) = (r*t)^n * exp(-r*t) / n!
static double factorial(int n) {
    // naive factorial for small n
    double f = 1.0;
    for(int i=1; i<=n; i++) f *= (double)i;
    return f;
}
double Psp(int n, double r, double t)
{
    // watch for large n. If n=0 or 1, fine. If you need bigger, handle carefully.
    double rt = r*t;
    double top = std::pow(rt, (double)n)*std::exp(-rt);
    double bot = factorial(n);
    return top/bot;
}

////////////////////////////////////////////////////
// 7. Probability of Detector "Click" Patterns
////////////////////////////////////////////////////

double Pm1010(double PG0, double PG1, double PG2,
              double PD0, double PD1, double PD2)
{
    return PG2*(PD0 + 2.0*PD1 + PD2) + 2.0*PG1*(PD1+PD2) + PG0*PD2;
}

double EtaTot(double PG0, double PG1, double PG2,
              double PD0, double PD1, double PD2)
{
    return 4.0 * Pm1010(PG0, PG1, PG2, PD0, PD1, PD2);
}

// Probability that a legitimate coincidence is present given the success signature
double PS(double PG0, double PG1, double PG2,
          double PD0, double PD1, double PD2)
{
    double numerator   = PG2*(PD0 + 2.0*PD1 + PD2);
    double denominator = Pm1010(PG0, PG1, PG2, PD0, PD1, PD2);
    if(denominator <= 0.0) return 0.0;
    return numerator/denominator;
}

// Overall Final Fidelity
double Fval(double PG0, double PG1, double PG2,
            double PD0, double PD1, double PD2,
            double x0, double sigma, double tmin, double tmax)
{
    double psval = PS(PG0, PG1, PG2, PD0, PD1, PD2);
    double fic   = Fic(x0, sigma, tmin, tmax);
    return psval*fic + (1.0 - psval)*(0.25);
}

////////////////////////////////////////////////////
// 8. Auxiliary: PG0, PG1, PG2 and PD0, PD1, PD2
////////////////////////////////////////////////////

// PG0[etaA, etaB] = (1 - etaA)*(1 - etaB)
double PG0(double etaA, double etaB)
{
    return (1.0 - etaA)*(1.0 - etaB);
}

double PG1(double etaA, double etaB)
{
    double firstTerm = (etaA*(1.0-etaB) + etaB*(1.0-etaA))*0.25;
    double secondTerm= 0.0625*etaA*etaB; // = 1/16
    return firstTerm + secondTerm;
}

double PG2(double etaA, double etaB)
{
    return 0.125*etaA*etaB; // 1/8
}

double PD0(double psp0)
{
    return std::pow(psp0, 4.0);
}

double PD1(double psp0)
{
    return std::pow(psp0,3.0)*(1.0 - psp0);
}

double PD2(double psp0)
{
    double t1 = std::pow(psp0,2.0);
    double t2 = std::pow(1.0-psp0,2.0);
    return t1*t2;
}

////////////////////////////////////////////////////
// 9. Putting it all together for day/night
////////////////////////////////////////////////////

double EtaTotDay(double RA, double w0, double h,
                 double theta, double theta2_,
                 double lambda, double sigma_tr, double eta_m,
                 double x0, double sigma, double tmin, double tmax,
                 double CD, double thetaFOV, double DeltaLambda)
{
    double eA = etaA(RA, w0, h, theta,  lambda, sigma_tr, eta_m, x0, sigma, tmin, tmax);
    double eB = etaB(RA, w0, h, theta2_,lambda, sigma_tr, eta_m, x0, sigma, tmin, tmax);
    double pg0 = PG0(eA, eB);
    double pg1 = PG1(eA, eB);
    double pg2 = PG2(eA, eB);

    double rd  = rday(CD, h, theta, eta_m, RA, thetaFOV, DeltaLambda);
    double dt  = (tmax - tmin);

    double p0 = Psp(0, rd, dt);
    double pd0= PD0(p0);
    double pd1= PD1(p0);
    double pd2= PD2(p0);

    return EtaTot(pg0, pg1, pg2, pd0, pd1, pd2);
}

double Fday(double RA, double w0, double h,
            double theta, double theta2_,
            double lambda, double sigma_tr, double eta_m,
            double x0, double sigma, double tmin, double tmax,
            double CD, double thetaFOV, double DeltaLambda)
{
    double eA = etaA(RA, w0, h, theta,  lambda, sigma_tr, eta_m, x0, sigma, tmin, tmax);
    double eB = etaB(RA, w0, h, theta2_,lambda, sigma_tr, eta_m, x0, sigma, tmin, tmax);
    double pg0 = PG0(eA, eB);
    double pg1 = PG1(eA, eB);
    double pg2 = PG2(eA, eB);

    double rd  = rday(CD, h, theta, eta_m, RA, thetaFOV, DeltaLambda);
    double dt  = (tmax - tmin);

    double p0  = Psp(0, rd, dt);
    double pd0 = PD0(p0);
    double pd1 = PD1(p0);
    double pd2 = PD2(p0);

    return Fval(pg0, pg1, pg2, pd0, pd1, pd2, x0, sigma, tmin, tmax);
}

// Night
double EtaTotNight(double RA, double w0, double h,
                   double theta, double theta2_,
                   double lambda, double sigma_tr, double eta_m,
                   double x0, double sigma, double tmin, double tmax,
                   double CD, double thetaFOV, double DeltaLambda,
                   double T)
{
    double eA = etaA(RA, w0, h, theta,  lambda, sigma_tr, eta_m, x0, sigma, tmin, tmax);
    double eB = etaB(RA, w0, h, theta2_,lambda, sigma_tr, eta_m, x0, sigma, tmin, tmax);
    double pg0 = PG0(eA, eB);
    double pg1 = PG1(eA, eB);
    double pg2 = PG2(eA, eB);

    double rn  = rnight(CD, h, theta, eta_m, RA, thetaFOV, DeltaLambda, lambda, T);
    double dt  = (tmax - tmin);

    double p0  = Psp(0, rn, dt);
    double pd0 = PD0(p0);
    double pd1 = PD1(p0);
    double pd2 = PD2(p0);

    return EtaTot(pg0, pg1, pg2, pd0, pd1, pd2);
}

double Fnight(double RA, double w0, double h,
              double theta, double theta2_,
              double lambda, double sigma_tr, double eta_m,
              double x0, double sigma, double tmin, double tmax,
              double CD, double thetaFOV, double DeltaLambda,
              double T)
{
    double eA = etaA(RA, w0, h, theta,  lambda, sigma_tr, eta_m, x0, sigma, tmin, tmax);
    double eB = etaB(RA, w0, h, theta2_,lambda, sigma_tr, eta_m, x0, sigma, tmin, tmax);
    double pg0 = PG0(eA, eB);
    double pg1 = PG1(eA, eB);
    double pg2 = PG2(eA, eB);

    double rn  = rnight(CD, h, theta, eta_m, RA, thetaFOV, DeltaLambda, lambda, T);
    double dt  = (tmax - tmin);

    double p0  = Psp(0, rn, dt);
    double pd0 = PD0(p0);
    double pd1 = PD1(p0);
    double pd2 = PD2(p0);

    return Fval(pg0, pg1, pg2, pd0, pd1, pd2, x0, sigma, tmin, tmax);
}


////////////////////////////////////////////////////
// Example main() showing how to use the routines
////////////////////////////////////////////////////
int main()
{
    // Just show a couple of test calls:

    // Example parameters:
    double CD    = 1500.0;
    double hVal  = 500e3;      // 500 km
    double DG    = 100e3;      // 100 km separation
    double RA    = 0.75;
    double w0    = 0.025;
    double thetaFOV = 1.0e-5;
    double DeltaLambda = 1.0e-9;
    double lambda_     = 8.0e-7;    // 800 nm
    double T_          = 300.0;     // K
    double cSpeed      = 3.0e8;
    double eta_m       = 0.25;
    double sigma_t     = 10e-9*cSpeed;  // wavepacket width in meters
    double tmin_       = -(40e-9)/2 + 10e-9;   // gating window
    double tmax_       =  (40e-9)/2 + 10e-9;
    double x0_         = cSpeed*(3e-9);  // path difference
    double sigma_tr    = 0.1;       // tracking error (m)

    // The "equidistant" angles:
    double th1 = theta_e(hVal, DG);
    double th2 = th1; // symmetrical

    // Evaluate night-time success probability:
    double etaNight = EtaTotNight(RA, w0, hVal, th1, th2,
                                  lambda_, sigma_tr, eta_m,
                                  x0_, sigma_t, tmin_, tmax_,
                                  CD, thetaFOV, DeltaLambda,
                                  T_);
    double fNight   = Fnight(RA, w0, hVal, th1, th2,
                             lambda_, sigma_tr, eta_m,
                             x0_, sigma_t, tmin_, tmax_,
                             CD, thetaFOV, DeltaLambda,
                             T_);

    std::cout << "Nighttime Success Probability = " << etaNight << "\n";
    std::cout << "Nighttime Fidelity           = " << fNight   << "\n";

    // Evaluate day-time success probability:
    double etaDay = EtaTotDay(RA, w0, hVal, th1, th2,
                              lambda_, sigma_tr, eta_m,
                              x0_, sigma_t, tmin_, tmax_,
                              CD, thetaFOV, DeltaLambda);
    double fDay  = Fday(RA, w0, hVal, th1, th2,
                        lambda_, sigma_tr, eta_m,
                        x0_, sigma_t, tmin_, tmax_,
                        CD, thetaFOV, DeltaLambda);

    std::cout << "Daytime  Success Probability = " << etaDay << "\n";
    std::cout << "Daytime  Fidelity           = " << fDay   << "\n";

    return 0;
}
