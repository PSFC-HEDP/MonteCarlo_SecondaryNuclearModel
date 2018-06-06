package MonteCarloParticleTracer;

import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;

public class Logger {

    private ArrayList<String> logs = new ArrayList<>();
    private boolean verbose = false;

    private HashMap<String, ArrayList<Long>> timers = new HashMap<>();

    public Logger() {
    }

    public Logger(boolean verbose) {
        this.verbose = verbose;
    }

    void addLog(String message){

        DateFormat dateFormat = new SimpleDateFormat("HH:mm:ss yyyy-MM-dd");
        message = dateFormat.format(new Date()) + ": " + message;

        logs.add(message);
        if (verbose) {
            System.out.println(message);
        }

    }

    void logTaskCompletion(){

        addLog("Done!");
        logs.add("");

    }

    void dumpLogToConsole(){
        System.out.println(toString());
    }

    void startTimer(String key){

        // Make this timer if it doesn't exists
        timers.putIfAbsent(key, new ArrayList<>());

        // Get the corresponding array list
        ArrayList<Long> times = timers.get(key);

        // Add the current time to the end of the array list
        times.add(System.nanoTime());

        // Put the array list back
        timers.put(key, times);

    }

    void stopTimer(String key){

        // Get the corresponding array list
        ArrayList<Long> times = timers.get(key);

        // The total time is now - the entry we stored when we started
        long totalTime = System.nanoTime() - times.get(times.size()-1);

        // Update the value to represent a total time
        times.set(times.size()-1, totalTime);

        // Put the array list back
        timers.put(key, times);

    }

    String getTimerResults(){

        StringBuilder results = new StringBuilder();
        results.append("Timer Name; Number Samples; Min Index; Max Index; Min; Average; Max; Total; Min/Average; Max/Average\n");

        // For each timer
        for (String key : timers.keySet()){

            ArrayList<Long> times = timers.get(key);
            double total = 0.0;
            double averageValue = 0.0;

            double maxValue = Long.MIN_VALUE;
            double minValue = Long.MAX_VALUE;

            int maxIndex = 0;
            int minIndex = 0;

            for (int i = 0; i < times.size(); i++){

                long time = times.get(i);
                total += times.get(i);

                if (maxValue < time) {
                    maxIndex = i;
                    maxValue = time;
                }

                if (minValue > time){
                    minIndex = i;
                    minValue = time;
                }

            }
            averageValue = total / times.size();


            results.append(key).append("; ");
            results.append(times.size()).append("; ");
            results.append(minIndex).append("; ");
            results.append(maxIndex).append("; ");
            results.append(formatTimeString(minValue));
            results.append(formatTimeString(averageValue));
            results.append(formatTimeString(maxValue));
            results.append(formatTimeString(total));
            results.append(minValue / averageValue).append(" ;");
            results.append(maxValue / averageValue).append(" ;\n");

        }

        return results.toString();
    }

    private String formatTimeString(double valueInNanoSeconds){

        String units;
        double conversion;

        if (valueInNanoSeconds > 1e9){
            units = "s";
            conversion = 1e-9;
        }
        else if (valueInNanoSeconds > 1e6){
            units = " ms";
            conversion = 1e-6;
        }
        else if (valueInNanoSeconds > 1e3){
            units = "us";
            conversion = 1e-3;
        }
        else{
            units = "ns";
            conversion = 1.0;
        }

        return String.format("%.2f %s; ", conversion*valueInNanoSeconds, units);

    }

    public String toString(){
        StringBuilder string = new StringBuilder();
        for (String log : logs){
            string.append(log).append("\n");
        }

        string.append("\n").append(getTimerResults());

        return string.toString();
    }

}
