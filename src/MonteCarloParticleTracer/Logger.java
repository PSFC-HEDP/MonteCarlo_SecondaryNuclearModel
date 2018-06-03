package MonteCarloParticleTracer;

import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;

public class Logger {

    private ArrayList<String> logs = new ArrayList<>();
    private boolean verbose = false;

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

        for (String log : logs) {
            System.out.println(log);
        }

    }

    public String toString(){
        String string = "";
        for (String log : logs){
            string += log + "\n";
        }
        return string;
    }

}
