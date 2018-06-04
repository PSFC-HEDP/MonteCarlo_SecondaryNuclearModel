package MonteCarloParticleTracer;

public class Exceptions {

    static class NoSourceInformationException extends Exception {
        NoSourceInformationException() {
            super("No source information was provided for this model");
        }
    }


    static class NoPlasmaSpecifiedException extends Exception {
        NoPlasmaSpecifiedException() {
            super("No plasma has been specified for this model");
        }
    }


    static class InvalidNumberProcessorsRequestedException extends Exception {
        InvalidNumberProcessorsRequestedException(int requested, int available){
            super(String.format("Invalid request for %d processors. This machine has %d available", requested, available));
        }
    }

    static class UnsupportedOperatingSystemException extends Exception {
        UnsupportedOperatingSystemException(String osName){
            super("No cStopPow libraries compiled for " + osName);
        }
    }


}
