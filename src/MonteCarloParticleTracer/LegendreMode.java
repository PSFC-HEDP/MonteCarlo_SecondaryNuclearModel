package MonteCarloParticleTracer;

import org.apache.commons.math3.util.CombinatoricsUtils;
import org.apache.commons.math3.util.FastMath;

public class LegendreMode{
    private int el, m;
    private double norm;
    private double magnitude;       // um

    LegendreMode(int el, int m, double magnitude) {
        this.el = el;
        this.m = m;
        this.magnitude = magnitude;

        this.norm = CombinatoricsUtils.factorial(el - m);
        norm /= CombinatoricsUtils.factorial(el + m);
        norm *= (2*el + 1);
        norm /= (4 * Math.PI);
        norm = Math.sqrt(norm);
    }

    double evaluate(double theta, double phi){
        double value = this.magnitude;
        value *= this.norm;
        value *= FastMath.cos(m * phi);
        value *= associatedLegendrePolyValue(theta);
        return value;
    }

    private double associatedLegendrePolyValue(double theta){
        double x = Math.cos(theta);     // Evaluation point
        double p0 = 1.0;
        int el = 0, m = 0;

        // Recurse until we find P(x, el=M, m=M)
        while ( m < this.m){
            p0 *= -(2*el + 1) * Math.sin(theta);
            m++;
            el++;
        }

        // If that's the one we need, return it
        if (el == this.el) return p0;

        // Else calculate P(x, el=M+1, m=M)
        double p1 = p0 * x * (2*el + 1);
        el++;

        // Recurse until we fine P(x, el=L, m=M)
        while ( el < this.el){
            double temp = p1;

            p1 = (2*el + 1)*x*p1 - (el + m)*p0;
            p1 /= (el - m + 1);

            p0 = temp;
            el++;
        }

        return p1;
    }

    double getMagnitude() {
        return magnitude;
    }

    int getEl() {
        return el;
    }

    int getM() {
        return m;
    }

    public double getNorm() {
        return norm;
    }

    void setMagnitude(double magnitude) {
        this.magnitude = magnitude;
    }


}