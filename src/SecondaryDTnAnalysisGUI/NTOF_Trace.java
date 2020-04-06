package SecondaryDTnAnalysisGUI;

import MonteCarloParticleTracer.*;
import PlottingAPI.Figure;
import PlottingAPI.LineProperties;

import java.io.File;

public class NTOF_Trace {

    double[] time;
    double[] voltage;

    public NTOF_Trace(File csvFile) {
        try {
            double[][] data = MonteCarloParticleTracer.Utils.parseCSV(csvFile);
            time    = data[0];
            voltage = data[1];
        }
        catch (Exception e){
            e.printStackTrace();
        }
    }

    public Figure plotData(double startTime, double endTime) {

        Figure figure = new Figure("", "Time (s)", "Voltage (V)");
        figure.plot(time, voltage, LineProperties.blackLine(2.0));
        figure.setXLimits(startTime, endTime);
        return figure;

    }
}
