//Homework 1

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <stdexcept>

using namespace std;

inline double sign(const double &a, const double &b) {
    return (b >= 0.0) ? fabs(a) : -fabs(a);
} 

ostream & operator<<(ostream & os, const vector<vector<double>> & vec);

class Atoms {
    private:
        string filename;

    public:
        Atoms(const string &filename);

        vector<vector<double>> read_xyz(const string &filename);

        friend ostream & operator<<(ostream & os, const vector<vector<double>> & vec);
};

class LennardJones : public Atoms {
    private:
        typedef vector<double> AtomCoord;
        typedef vector<vector<double>> Coordinates;

    public:
        double sigma = 2.951;
        double epsilon = 5.29;

        LennardJones(const string &filename);

        friend ostream & operator<<(ostream & os, const vector<vector<double>> & vec);

        double calculate_LJ(double r_ij);
        
        double calculate_distance(AtomCoord coord1, AtomCoord coord2);

        double calculate_total_energy(Coordinates coordinates);

        double calculate_pair_energy(Coordinates coordinates, int i_particle);

        void run_LJ(string &filename);
};

class FiniteDifference : public LennardJones {
    private:
        vector<double> h = {0.1, 0.01, 0.001, 0.0001};
        typedef vector<double> AtomCoord;
        typedef vector<vector<double>> Coordinates;

    public:
        FiniteDifference(const string &filename);

        friend ostream & operator<<(ostream & os, Coordinates vec);

        void printFormattedData(const string & label, Coordinates data, double stepsize);

        vector<vector<double>> analytical_force(Coordinates coordinates);

        double forward_difference(Coordinates coordinates, int i_particle, int dim, double h);

        double central_difference(Coordinates coordinates, int i_particle, int dim, double h);

        void steepest_descent(Coordinates coordinates);

        void run_FD(string &filename);
};

struct Bracketmethod {
    double ax, bx, cx, fa, fb, fc;
    template <class T>
    void bracket(const double a, const double b, T &func) {
        const double GOLD=1.618034, GLIMIT=100.0, TINY=1.0e-20;
        ax=a; bx=b;
        double fu;
        fa=func(ax);
        fb=func(bx);
        if (fb > fa) {
            swap(ax, bx);
            swap(fb, fa);
        }
        cx=bx+GOLD*(bx-ax);
        fc=func(cx);
        while (fb > fc) {
            double r=(bx-ax)*(fb-fc);
            double q=(bx-cx)*(fb-fa);
            double u=bx-((bx-cx)*q-(bx-ax)*r)/(2.0*sign(max(abs(q-r),TINY),q-r));
            double ulim=bx+GLIMIT*(cx-bx);
            if ((bx-u)*(u-cx) > 0.0) {
                fu=func(u);
                if (fu < fc) {
                    ax=bx;
                    bx=u;
                    fa=fb;
                    fb=fu;
                    return; 
                } 
                else if (fu > fb) {
                    cx=u;
                    fc=fu;
                    return;
                }
                u=cx+GOLD*(cx-bx);
                fu=func(u);
            } 
            else if ((cx-u)*(u-ulim) > 0.0) {
                fu=func(u);
                if (fu < fc) {
                    shft3(bx,cx,u,u+GOLD*(u-cx));
                    shft3(fb,fc,fu,func(u));
                }
            } 
            else if ((u-ulim)*(ulim-cx) >= 0.0) {
                u=ulim; 
                fu=func(u);
            } 
            else {
                u=cx+GOLD*(cx-bx);
                fu=func(u);
            }
            shft3(ax, bx, cx, u);
            shft3(fa, fb, fc, fu);
        } 
    }
    inline void shft2(double &a, double &b, const double c) {
        a=b;
        b=c;
    }
    inline void shft3(double &a, double &b, double &c, const double d) {
        a=b;
        b=c;
        c=d;
    }
    inline void mov3(double &a, double &b, double &c, const double d, const double e, const double f) {
        a=d;
        b=e;
        c=f; 
    }
};

struct Golden : public Bracketmethod {
    double xmin, fmin;
    const double tol;
    Golden(const double toll=3.0e-8) : tol(toll) {}
    template <class T>
    double minimize(T &func) {
        const double R=0.61803399,C=1.0-R;
        double x1, x2;
        double x0=ax;
        double x3=cx;
        if (abs(cx-bx) > abs(bx-ax)) {
            x1=bx;
            x2=bx+C*(cx-bx);
        }  
        else {
            x2=bx;
            x1=bx-C*(bx-ax);
        }
        double f1=func(x1);
        double f2=func(x2);
        while (abs(x3-x0) > tol*(abs(x1)+abs(x2))) {
            if (f2 < f1) {
            shft3(x0, x1, x2, R*x2+C*x3); 
            shft2(f1, f2, func(x2));
            }
            else {
            shft3(x3, x2, x1, R*x1+C*x0);
            shft2(f2, f1, func(x1));
            }
        }
        if (f1 < f2) {
            xmin=x1;
            fmin=f1;
        } 
        else {
            xmin=x2;
            fmin=f2;
        }
        return xmin;
    }
};

/*
struct Brent : public Bracketmethod {
    double xmin, fmin;
    const double tol;
    Brent(const double toll=3.0e-8) : tol(toll) {}
    template <class T>
    double minimize(T &func) {
        const int ITMAX=100;
        const double CGOLD=0.3819660;
        const double ZEPS=numeric_limits<double>::epsilon()*1.0e-3;
        double a, b, d=0.0, etemp, fu, fv, fw, fx;
        double p, q, r, tol1, tol2, u, v, w, x, xm;
        double e=0.0;
        a=(ax<cx ? ax:cx);
        b=(ax>cx ? ax:cx);
        x=w=v=bx;
        fw=fv=fx=func(x);
        for (int iter=0; iter<ITMAX; iter++) {
            xm=0.5*(a+b);
            tol2=2.0*(tol1=tol*abs(x)+ZEPS);
            if (abs(x-xm) <= (tol2-0.5*(b-a))) { 
                fmin=fx;
                return xmin=x;
            }
            if (abs(e) > tol1) {
                r=(x-w)*(fx-fv);
                q=(x-v)*(fx-fw);
                p=(x-v)*q-(x-w)*r;
                q=2.0*(q-r);
                if (q > 0.0) p =-p;
                q=abs(q);
                etemp=e;
                e=d;
                if (abs(p) >= abs(0.5*q*etemp) || p <= q*(a-x)
                || p >= q*(b-x))
                d=CGOLD*(e=(x >= xm ? a-x : b-x));
                else {
                    d=p/q;
                    u=x+d;
                    if (u-a < tol2 || b-u < tol2)
                    d=SIGN(tol1,xm-x);
                }
            } 
            else {
                d=CGOLD*(e=(x >= xm ? a-x : b-x));
            }
            u=(abs(d) >= tol1 ? x+d : x+sign(tol1,d));
            fu=func(u);
            if (fu <= fx) {
                if (u >= x) a=x; else b=x;
                shft3(v,w,x,u);
                shft3(fv,fw,fx,fu);
            } 
            else {
                if (u < x) a=u; else b=u;
                if (fu <= fw || w == x) {
                    v=w;
                    w=u;
                    fv=fw;
                    fw=fu;
                } 
                else if (fu <= fv || v == x || v == w) {
                    v=u;
                    fv=fu;
                }
            }
        }
        throw ("Too many iterations in brent");
    }
}; */
