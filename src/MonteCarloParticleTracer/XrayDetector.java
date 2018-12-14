package MonteCarloParticleTracer;

import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.exception.OutOfRangeException;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.util.FastMath;

import java.io.File;
import java.util.ArrayList;
import java.util.Scanner;

public class XrayDetector {

    Vector3D pinholeCenter;           // Physical location of the pinhole (cm)
    double pinholeDiameter;             // Diameter of the pinhole (cm)

    Vector3D detectorCenter;          // Physical location of the detector's center (cm)
    Vector3D detectorBasisBasisVectorX; // "X" unit basis vector in the detector plane
    Vector3D detectorPlaneBasisVectorY; // "Y" unit basis vector in the detector plane

    // Array list of all of the filters in front of the detector
    ArrayList<Filter> filters = new ArrayList<>();


    public XrayDetector(Vector3D pinholeCenter, Vector3D detectorCenter, double pinholeDiameter) {

        // Set the fields
        this.pinholeCenter = pinholeCenter;
        this.detectorCenter = detectorCenter;
        this.pinholeDiameter = pinholeDiameter;

        // Generate a basis set for the detector plane
        Vector3D norm = detectorCenter.normalize();
        detectorBasisBasisVectorX = norm.orthogonal();
        detectorPlaneBasisVectorY = Vector3D.crossProduct(norm, detectorBasisBasisVectorX);

    }


    public void addFilter(double density, double thickness, File massAttenuationFile){
        filters.add(new Filter(density, thickness, massAttenuationFile));
    }

    public void addGraphiteFilter(double thickness){
        addFilter(2.266, thickness, DataFiles.C_NIST_MassAttenuationFile);
    }

    public void addKaptonFilter(double thickness){
        addFilter(1.43, thickness, DataFiles.Kapton_NIST_MassAttenuationFile);
    }

    public void addAlFilter(double thickness){
        addFilter(2.7, thickness, DataFiles.Al_NIST_MassAttenuationFile);
    }

    public void addCuFilter(double thickness){
        addFilter(8.96, thickness, DataFiles.Cu_NIST_MassAttenuationFile);
    }

    public void addVFilter(double thickness){
        addFilter(5.8, thickness, DataFiles.V_NIST_MassAttenuationFile);
    }

    public void addGeFilter(double thickness){
        addFilter(5.323, thickness, DataFiles.Ge_NIST_MassAttenuationFile);
    }

    public void addMoFilter(double thickness){
        addFilter(10.28, thickness, DataFiles.Mo_NIST_MassAttenuationFile);
    }


    public double getTransmission(double energy){
        double transmission = 1.0;
        for (Filter filter : filters){
            transmission *= filter.getTransmission(energy);
        }
        return transmission;
    }


    public Tally_2D generateImage(double[] xBins, double[] yBins, Plasma plasma){

        final int NUM_X_RAYS = (int) 4e7;


        // Init the tally
        Tally_2D image = new Tally_2D(xBins, yBins);

        // Get the radial distribution
        Distribution radialSourceDistribution = plasma.getRadialBremDistribution();

        // Get the polar distribution from the plasma before we start
        Distribution polarSourceDistribution = null;
        if (plasma.boundsVaryInTheta()){
            polarSourceDistribution = plasma.getPolarDistribution();
        }


        // Loop through all of the x-rays we want to sample
        for (int i = 0; i < NUM_X_RAYS; i++){

            // ********************************
            // Sample the x-ray source position
            // ********************************

            // Start by sampling a spherically random normal vector
            Vector3D sourcePosition = Utils.sampleRandomNormalizedVector();

            // Convert it to spherical
            double[] coordinates = Utils.getSphericalFromVector(sourcePosition);
            double theta = coordinates[1];
            double phi   = coordinates[2];

            // If we have a polar distribution, use it instead
            if (polarSourceDistribution != null){
                theta = polarSourceDistribution.sample();
                sourcePosition = Utils.getVectorFromSpherical(1.0, theta, phi);
            }

            // Use our radial distribution and the plasma bounds to determine the magnitude of our position vector
            double rNorm = radialSourceDistribution.sample();
            double rMax  = plasma.getOuterRadiusBound(theta, phi);
            double rMin  = plasma.getInnerRadiusBound(theta, phi);

            double r = (rMax - rMin) * rNorm + rMin;

            // Rescale our position vector to this magnitude
            sourcePosition = sourcePosition.scalarMultiply(r);

            // Get the electron temperature at this position
            double Te = plasma.getElectronTemperature(sourcePosition);


            // ******************************************************
            // Sample where the x-ray passes through the pinhole area
            // ******************************************************

            double rPinhole = (pinholeDiameter / 2.0) * Math.sqrt(FastMath.random());
            double thetaPinhole = 2 * Math.PI * FastMath.random();

            // Create the vector that points to where the x-ray is passing through the pinhole
            // r_location = r_center + R * ( cos(theta) e_x + sin(theta) e_y )

            Vector3D pinholeLocation = detectorBasisBasisVectorX.scalarMultiply(FastMath.cos(thetaPinhole));
            pinholeLocation = pinholeLocation.add(detectorPlaneBasisVectorY.scalarMultiply(FastMath.sin(thetaPinhole)));
            pinholeLocation = pinholeLocation.scalarMultiply(rPinhole);
            pinholeLocation = pinholeLocation.add(pinholeCenter);


            // *******************************************
            // Determine where the x-ray hits the detector
            // *******************************************

            // Create the vector that points from the start to the pinhole
            Vector3D pathVector = pinholeLocation.subtract(sourcePosition);

            // Scale this vector to the detector plane
            double factor = Math.pow(detectorCenter.getNorm(), 2);
            factor /= pathVector.dotProduct(detectorCenter);
            pathVector = pathVector.scalarMultiply(factor);

            // Calculate the "x, y" of this vector
            double x = pathVector.dotProduct(detectorBasisBasisVectorX);
            double y = pathVector.dotProduct(detectorPlaneBasisVectorY);


            // ***********************
            // Sample the x-ray energy
            // ***********************

            // This samples the exp according to the local Te
            double energy = -Te * Math.log(1 - FastMath.random());


            // *********************
            // Add this to the tally
            // *********************

            double weight = getTransmission(energy);
            image.addValue(x, y, weight);

        }

        return image;

    }

    public double[][] generateImage_Old(double[] xBins, double[] yBins, Plasma plasma) throws Exception{

        // ****************************************
        // Constants for the plasma volume integral
        // ****************************************

        final int NUM_PHI_NODES    = 50;
        final int NUM_THETA_NODES  = 50;
        final int NUM_RADIUS_NODES = 50;

        final double MIN_PHI = 0.0;
        final double MAX_PHI = 2 * Math.PI;

        final double MIN_THETA = 0.0;
        final double MAX_THETA = 1 * Math.PI;

        double[][] results = new double[xBins.length-1][yBins.length-1];


        // *********************************
        // Loop over x image pixel positions
        // *********************************

        for (int xIndex = 1; xIndex < xBins.length; xIndex++){

            double xPosMagnitude = 0.5*(xBins[xIndex] + xBins[xIndex - 1]);


            // *********************************
            // Loop over y image pixel positions
            // *********************************

            for (int yIndex = 1; yIndex < yBins.length; yIndex++){

                double yPosMagnitude = 0.5*(yBins[yIndex] + yBins[yIndex - 1]);

                // Make a vector to the pixel
                Vector3D pixelLocation = new Vector3D(0.0, 0.0, 0.0);
                pixelLocation = pixelLocation.add(detectorCenter);
                pixelLocation = pixelLocation.add(detectorBasisBasisVectorX.scalarMultiply(xPosMagnitude));
                pixelLocation = pixelLocation.add(detectorPlaneBasisVectorY.scalarMultiply(yPosMagnitude));

                //System.out.printf("%d %d\n", xIndex, yIndex);


                // ************************************************
                // Loop over the phi dimension of the source plasma
                // ************************************************

                for (int phiIndex = 0; phiIndex < NUM_PHI_NODES; phiIndex++) {

                    double phi_1 = MIN_PHI + (phiIndex + 1) * (MAX_PHI - MIN_PHI) / (NUM_PHI_NODES);
                    double phi_0 = MIN_PHI + (phiIndex + 0) * (MAX_PHI - MIN_PHI) / (NUM_PHI_NODES);
                    double dPhi = (phi_1 - phi_0);


                    // **************************************************
                    // Loop over the theta dimension of the source plasma
                    // **************************************************

                    for (int thetaIndex = 0; thetaIndex < NUM_THETA_NODES; thetaIndex++) {

                        double theta_1 = MIN_THETA + (thetaIndex + 1) * (MAX_THETA - MIN_THETA) / (NUM_THETA_NODES);
                        double theta_0 = MIN_THETA + (thetaIndex + 0) * (MAX_THETA - MIN_THETA) / (NUM_THETA_NODES);
                        double dTheta = (theta_1 - theta_0);

                        // Evaluate all of the outer radius bounds corresponding to this solid angle slice
                        double rMax_11 = plasma.getOuterRadiusBound(theta_1, phi_1);
                        double rMax_10 = plasma.getOuterRadiusBound(theta_1, phi_0);
                        double rMax_01 = plasma.getOuterRadiusBound(theta_0, phi_1);
                        double rMax_00 = plasma.getOuterRadiusBound(theta_0, phi_0);

                        // Evaluate all of the inner radius bounds corresponding to this solid angle slice
                        double rMin_11 = plasma.getInnerRadiusBound(theta_1, phi_1);
                        double rMin_10 = plasma.getInnerRadiusBound(theta_1, phi_0);
                        double rMin_01 = plasma.getInnerRadiusBound(theta_0, phi_1);
                        double rMin_00 = plasma.getInnerRadiusBound(theta_0, phi_0);


                        // ***************************************************
                        // Loop over the radius dimension of the source plasma
                        // ***************************************************

                        for (int radiusIndex = 0; radiusIndex < NUM_RADIUS_NODES; radiusIndex++) {

                            // Determine the radius binEdges corresponding to theta_1, phi_1
                            double r_111 = rMin_11 + (radiusIndex + 1) * (rMax_11 - rMin_11) / NUM_RADIUS_NODES;
                            double r_011 = rMin_11 + (radiusIndex + 0) * (rMax_11 - rMin_11) / NUM_RADIUS_NODES;
                            double dR_11 = (r_111 - r_011);

                            // Determine the radius binEdges corresponding to theta_1, phi_0
                            double r_110 = rMin_10 + (radiusIndex + 1) * (rMax_10 - rMin_10) / NUM_RADIUS_NODES;
                            double r_010 = rMin_10 + (radiusIndex + 0) * (rMax_10 - rMin_10) / NUM_RADIUS_NODES;
                            double dR_10 = (r_110 - r_010);

                            // Determine the radius binEdges corresponding to theta_0, phi_1
                            double r_101 = rMin_01 + (radiusIndex + 1) * (rMax_01 - rMin_01) / NUM_RADIUS_NODES;
                            double r_001 = rMin_01 + (radiusIndex + 0) * (rMax_01 - rMin_01) / NUM_RADIUS_NODES;
                            double dR_01 = (r_101 - r_001);

                            // Determine the radius binEdges corresponding to theta_0, phi_0
                            double r_100 = rMin_00 + (radiusIndex + 1) * (rMax_00 - rMin_00) / NUM_RADIUS_NODES;
                            double r_000 = rMin_00 + (radiusIndex + 0) * (rMax_00 - rMin_00) / NUM_RADIUS_NODES;
                            double dR_00 = (r_100 - r_000);


                            // Create a vector for all 8 binEdges
                            Vector3D vec_111 = Utils.getVectorFromSpherical(r_111, theta_1, phi_1);
                            Vector3D vec_011 = Utils.getVectorFromSpherical(r_011, theta_1, phi_1);

                            Vector3D vec_110 = Utils.getVectorFromSpherical(r_110, theta_1, phi_0);
                            Vector3D vec_010 = Utils.getVectorFromSpherical(r_010, theta_1, phi_0);

                            Vector3D vec_101 = Utils.getVectorFromSpherical(r_101, theta_0, phi_1);
                            Vector3D vec_001 = Utils.getVectorFromSpherical(r_001, theta_0, phi_1);

                            Vector3D vec_100 = Utils.getVectorFromSpherical(r_100, theta_0, phi_0);
                            Vector3D vec_000 = Utils.getVectorFromSpherical(r_000, theta_0, phi_0);


                            // Throw all the node quantities in arrays so we can loop through them easily
                            double[] rNodes = new double[]{r_111, r_011, r_110, r_010, r_101, r_001, r_100, r_000};
                            double[] dRs = new double[]{dR_11, dR_11, dR_10, dR_10, dR_01, dR_01, dR_00, dR_00};
                            double[] thetaNodes = new double[]{theta_1, theta_1, theta_1, theta_1, theta_0, theta_0, theta_0, theta_0};
                            Vector3D[] vectors = new Vector3D[]{vec_111, vec_011, vec_110, vec_010, vec_101, vec_001, vec_100, vec_000};


                            Vector3D sourceVector = new Vector3D(0.0, 0.0, 0.0);          // Vector in the center of the volume element
                            double elementVolume = 0.0;                                          // Volume of the element


                            // **************************************
                            // Loop over the 8 resultant source nodes
                            // **************************************

                            for (int nodeIndex = 0; nodeIndex < rNodes.length; nodeIndex++) {

                                sourceVector = sourceVector.add(vectors[nodeIndex].scalarMultiply(1.0 / 8.0));

                                double value = rNodes[nodeIndex] * rNodes[nodeIndex];       // r^2 from Jacobian
                                value *= FastMath.sin(thetaNodes[nodeIndex]);               // sin(theta) from Jacobian
                                value *= dRs[nodeIndex] * dTheta * dPhi / 8.0;              // Integration widths

                                elementVolume += value;

                            }

                            // Calculate the vector between the source vector and the pixel vector
                            Vector3D pathVector = pixelLocation.subtract(sourceVector);


                            // Scale the path vector to be in the plane of the pinhole
                            double factor = pinholeCenter.getNorm();
                            factor -= sourceVector.dotProduct(pinholeCenter.normalize());
                            factor /= pathVector.dotProduct(pinholeCenter.normalize());

                            // Find the vector that points to where this path intersects the pinhole plane
                            Vector3D pinholeIntersect = pathVector.scalarMultiply(factor);
                            pinholeIntersect = pinholeIntersect.add(sourceVector);

                            // Determine how far that intersect is away from the pinhole
                            double distance = pinholeIntersect.subtract(pinholeCenter).getNorm();

                            // If that distance is smaller than the radius, we pass through!
                            if (2.0*distance <= pinholeDiameter) {

                                // TODO: Calculate the x-ray emission from this source volume and add it to the pixel

                                double emission = Math.pow(elementVolume, 1);
                                results[xIndex-1][yIndex-1] += emission;


                            }
                        }
                    }
                }

                System.out.printf("%.4e, ", results[xIndex-1][yIndex-1]);

            }

            System.out.println();

        }

        return results;
    }

    class Filter {

        double density;             // g / cc
        double thickness;           // cm

        PolynomialSplineFunction logMassAttenFunction;

        public Filter(double density, double thickness, File massAttenuationFile) {
            this.density = density;
            this.thickness = thickness;

            try {

                Scanner s = new Scanner(massAttenuationFile);
                ArrayList<Double> energyList = new ArrayList<>();
                ArrayList<Double> coefficientList = new ArrayList<>();

                while (s.hasNextDouble()) {
                    energyList.add(s.nextDouble());         // MeV
                    coefficientList.add(s.nextDouble());    // cm^2 / g
                    s.nextDouble();
                }


                double[] energyLogs = new double[energyList.size()];
                for (int i = 0; i < energyLogs.length; i++){
                    energyLogs[i] = Math.log(energyList.get(i));
                }

                double[] coefficientLogs = new double[coefficientList.size()];
                for (int i = 0; i < coefficientLogs.length; i++){
                    coefficientLogs[i] = Math.log(coefficientList.get(i));
                }

                logMassAttenFunction = new LinearInterpolator().interpolate(energyLogs, coefficientLogs);

            }
            catch (Exception e){
                e.printStackTrace();
                System.exit(-1);
            }
        }

        public double getTransmission(double energy) {
            try {
                double massAtten = logMassAttenFunction.value(Math.log(0.001 * energy));
                massAtten = Math.exp(massAtten);
                return Math.exp(-massAtten * density * thickness);
            }
            catch (OutOfRangeException e){
                if(e.getArgument().doubleValue() > e.getHi().doubleValue()){
                    return 1.0;
                }
                else if (e.getArgument().doubleValue() < e.getLo().doubleValue()){
                    return  0.0;
                }

                return Double.NaN;
            }
        }
    }

}

