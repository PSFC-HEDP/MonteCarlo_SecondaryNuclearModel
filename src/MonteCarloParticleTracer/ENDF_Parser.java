package MonteCarloParticleTracer;

import com.sun.org.apache.bcel.internal.generic.IF_ACMPEQ;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

import java.io.File;
import java.util.ArrayList;
import java.util.Scanner;

public class ENDF_Parser {

    private final double NUETRON_MASS_AMU = 1.008664;
    private File ENDF_File;


    public ENDF_Parser(File ENDF_File) {
        this.ENDF_File = ENDF_File;
    }


    public void parseENDF(int requestedReaction) throws Exception{

        // Things that we might wanna return somehow
        int Z, A;
        double mass;
        double Q;

        PolynomialSplineFunction crossSection;
        ArrayList<PolynomialSplineFunction> legendreCoefficients = new ArrayList<>();


        // Open the file
        Scanner s = new Scanner(ENDF_File);


        // Loop the file
        while (s.hasNextLine()){

            // Get the entries
            String[] entries = getEntries(s.nextLine());
            int MAT  = Integer.parseInt(entries[6]);
            int MF   = Integer.parseInt(entries[7]);
            int MT   = Integer.parseInt(entries[8]);
            int LINE = Integer.parseInt(entries[9]);

            // We don't care about anything at isn't our reaction
            if (MT != requestedReaction){
                continue;
            }

            // XS Data
            if (MF == 3){

                // *********************************************
                // First line
                // ZAID (1000Z + A), AWR (mT/mn), 0, 0, 0, 0, ...
                // *********************************************

                Z = (int) Math.floor(toDouble(entries[0]) / 1000);
                A = (int) Math.round(toDouble(entries[0])) - 1000*Z;
                mass = NUETRON_MASS_AMU * toDouble(entries[1]);



                // *************************************************************************************
                // Second line
                // QM (Q value), QI (don't care), 0, LR (complex flag), NR (num regions), NP (num pairs)
                // *************************************************************************************

                entries = getEntries(s.nextLine());

                Q              = toDouble(entries[0]);
                double QI      = toDouble(entries[1]);
                int LR         = Integer.valueOf(entries[3]);
                int numRegions = Integer.valueOf(entries[4]);
                int numPoints  = Integer.valueOf(entries[5]);

                if (LR != 0)         System.err.println("LR == 1 not supported!");
                if (numRegions != 1) System.err.println("Multiple interpolation regions not supported!");



                // *************************************************************************************
                // Third line (This breaks if NR != 1)
                // NP(i) Interp(i)
                // *************************************************************************************

                entries = getEntries(s.nextLine());

                int interpIndex = Integer.valueOf(entries[1]);
                if (interpIndex != 2)   System.err.println("Only linear interpolation is supported!");



                // *******
                // XS Data
                // *******

                // Init the data
                double[] energies = new double[numPoints];      // MeV
                double[] values   = new double[numPoints];      // cm^2

                // Loop through the lines we need (there are 3 point pairs per line)
                int numLines = numPoints / 3;
                for (int i = 0; i < numLines; i++){
                    entries = getEntries(s.nextLine());

                    energies[3*i + 0] = toDouble(entries[0]) * 1e-6;
                    energies[3*i + 1] = toDouble(entries[2]) * 1e-6;
                    energies[3*i + 2] = toDouble(entries[4]) * 1e-6;

                    values  [3*i + 0] = toDouble(entries[1]) * 1e-24;
                    values  [3*i + 1] = toDouble(entries[3]) * 1e-24;
                    values  [3*i + 2] = toDouble(entries[5]) * 1e-24;
                }

                // We need to pick any remaining points on the following line
                int remainingPoints = Math.floorMod(numPoints, 3);
                if (remainingPoints != 0){

                    entries = getEntries(s.nextLine());
                    for (int i = 0; i < remainingPoints; i++){
                        energies[3*numLines + i] = toDouble(entries[2*i + 0]) * 1e-6;
                        values  [3*numLines + i] = toDouble(entries[2*i + 1]) * 1e-24;
                    }
                }

                // Build the cross section function
                crossSection = new LinearInterpolator().interpolate(energies, values);



            }


            // Angular data
            if (MF == 4) {

                // **************************************************************
                // First line
                // ZAID (1000Z + A), AWR (mT/mn), 0, LTT (format flag), 0, 0, ...
                // **************************************************************

                Z = (int) Math.floor(toDouble(entries[0]) / 1000);
                A = (int) Math.round(toDouble(entries[0])) - 1000*Z;
                mass = NUETRON_MASS_AMU * toDouble(entries[1]);
                int LTT = Integer.valueOf(entries[3]);

                if (LTT == 2)         System.err.println("Tabulated angle distributions not supported!");



                // **************************************************************
                // Second line
                // 0, AWR (mT/mn), LI (iso flag), LCT (Frame flag), 0, 0, ...
                // **************************************************************

                entries = getEntries(s.nextLine());

                mass = NUETRON_MASS_AMU * toDouble(entries[1]);
                int LI  = Integer.valueOf(entries[2]);
                int LCT = Integer.valueOf(entries[3]);

                if (LI == 1)         System.err.println("Isotropic reactions not supported!");
                if (LCT == 1)        System.err.println("Lab frame angular distributions not supported!");



                // **************************************************************
                // Third line
                // 0, 0, 0, 0, NR (num regions), NE (num energies), ...
                // **************************************************************

                entries = getEntries(s.nextLine());

                int numRegions  = Integer.valueOf(entries[4]);
                int numEnergies = Integer.valueOf(entries[5]);


                if (numRegions != 1) System.err.println("Multiple interpolation regions not supported!");



                // **************************************************************
                // Fourth line (This breaks if NR != 1)
                // NE(i), Interp(1)
                // **************************************************************

                entries = getEntries(s.nextLine());

                int interpIndex = Integer.valueOf(entries[1]);
                if (interpIndex != 2)   System.err.println("Only linear interpolation is supported!");



                // ****************
                // Coefficient data
                // ****************

                ArrayList<double[]> coefficients = new ArrayList<>();
                double[] energies = new double[numEnergies];        // MeV

                // Loop through the energies
                for (int i = 0; i < numEnergies; i++){

                    //  0, energy[i], 0, 0, NL (num coefficients), 0
                    entries = getEntries(s.nextLine());

                    energies[i] = toDouble(entries[1]) * 1e-6;
                    int numCoefficients = Integer.valueOf(entries[4]);


                    // Add more arrays to the ArrayList if needed
                    while (numCoefficients > coefficients.size()){
                        coefficients.add(new double[numEnergies]);
                    }


                    // Loop through the lines we need (there are up to 6 coefficients per line)
                    int numLines = numCoefficients / 6;
                    for (int j = 0; j < numLines; j++){
                        entries = getEntries(s.nextLine());

                        // Loop through the 6 entries
                        for (int k = 0; k < 6; k++){
                            double[] values = coefficients.get(6*j + k);
                            values[i] = toDouble(entries[k]);
                            coefficients.set(6*j + k, values);
                        }
                    }

                    // We need to pick any remaining points on the following line
                    int remainingCoefficients = Math.floorMod(numCoefficients, 6);
                    if (remainingCoefficients != 0){

                        entries = getEntries(s.nextLine());
                        for (int k = 0; k < remainingCoefficients; k++){

                            double[] values = coefficients.get(6*numLines + k);
                            values[i] = toDouble(entries[k]);
                            coefficients.set(6*numLines + k, values);
                        }
                    }
                }

                // Build the coefficient functions
                for (int i = 0; i < coefficients.size(); i++){
                    legendreCoefficients.add(new LinearInterpolator().interpolate(energies, coefficients.get(i)));
                }
            }
        }

        // Close the file
        s.close();

    }

    private String[] getEntries(String line){

        String[] entries = new String[10];

        // First 6 entries have the same length (of 11)
        for (int i = 0; i < 6; i++){
            entries[i] = line.substring(i*11, (i+1)*11).trim();
        }

        entries[6] = line.substring(67-1, 67-1 + 4).trim();        // MAT
        entries[7] = line.substring(71-1, 71-1 + 2).trim();        // MF
        entries[8] = line.substring(73-1, 73-1 + 3).trim();        // MT
        entries[9] = line.substring(76-1, 76-1 + 5).trim();        // Line Number

        return entries;
    }

    private Double toDouble(String entry){

        // Make sure this is in the format we expect
        if (!entry.contains("+") && !entry.contains("-")){
            System.err.println("Invalid Double Format!");
        }

        // Check the sign of the power and base
        boolean negativePower = !entry.contains("+");
        boolean negativeBase  = entry.startsWith("-");

        // Remove the "-" if needed
        if (negativeBase)   entry = entry.substring(1, entry.length());

        // Split the string into its two parts
        String[] temp = entry.split("\\+|-");

        // Get the base
        double base = Double.valueOf(temp[0]);
        if (negativeBase)   base *= -1;

        // Get the power
        int power   = Integer.valueOf(temp[1]);
        if (negativePower)  power *= -1;

        // Return
        return base * Math.pow(10, power);
    }

}
