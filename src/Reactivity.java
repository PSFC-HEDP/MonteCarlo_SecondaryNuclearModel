import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;

/**
 * Created by lahmann on 2016-06-16.
 */
public class Reactivity {

    final static Reactivity DDn_Reactivity = ddn();
    final static Reactivity DDp_Reactivity = ddp();
    final static Reactivity D3Hep_Reactivity = d3Hep();


    private PolynomialSplineFunction value;

    private static Reactivity ddn(){
        return new Reactivity(Utils.DDn_REACTIVITY_FILE);
    }

    private static Reactivity ddp(){
        return new Reactivity(Utils.DDp_REACTIVITY_FILE);
    }

    private static Reactivity d3Hep(){
        return new Reactivity(Utils.D3Hep_REACTIVITY_FILE);
    }

    public Reactivity(String dataFilename){
        this(new File(dataFilename));
    }

    public Reactivity(File dataFile){

        ArrayList<Double> T = new ArrayList<>();
        ArrayList<Double> sigV = new ArrayList<>();
        try{
            Scanner s = new Scanner(dataFile);

            while (s.hasNext()){
                T.add(s.nextDouble());
                sigV.add(s.nextDouble());
            }

            double[] eArray = new double[T.size()];
            double[] sigArray = new double[sigV.size()];
            for (int i = 0; i < T.size(); i++){
                eArray[i] = T.get(i);
                sigArray[i] = sigV.get(i);
            }

            value = new LinearInterpolator().interpolate(
                    eArray, sigArray);
        }
        catch (IOException e){
            e.printStackTrace();
        }
    }


    double evaluate(double temperature){
        return value.value(temperature);
    }
}
