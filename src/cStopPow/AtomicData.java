/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 3.0.10
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

package cStopPow;

public class AtomicData {
  private transient long swigCPtr;
  protected transient boolean swigCMemOwn;

  protected AtomicData(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  protected static long getCPtr(AtomicData obj) {
    return (obj == null) ? 0 : obj.swigCPtr;
  }

  protected void finalize() {
    delete();
  }

  public synchronized void delete() {
    if (swigCPtr != 0) {
      if (swigCMemOwn) {
        swigCMemOwn = false;
        cStopPowJNI.delete_AtomicData(swigCPtr);
      }
      swigCPtr = 0;
    }
  }

  public static double get_AMU(int Z) {
    return cStopPowJNI.AtomicData_get_AMU(Z);
  }

  public static double get_rho(int Z) {
    return cStopPowJNI.AtomicData_get_rho(Z);
  }

  public static String get_symbol(int Z) {
    return cStopPowJNI.AtomicData_get_symbol(Z);
  }

  public static int get_num_from_symbol(String symbol) {
    return cStopPowJNI.AtomicData_get_num_from_symbol(symbol);
  }

  public static String get_name(int Z) {
    return cStopPowJNI.AtomicData_get_name(Z);
  }

  public static int get_num_from_name(String name) {
    return cStopPowJNI.AtomicData_get_num_from_name(name);
  }

  public static double get_mean_ionization(int Z) {
    return cStopPowJNI.AtomicData_get_mean_ionization(Z);
  }

  public static SWIGTYPE_p_std__arrayT_double_5_t get_shell_coeff(int Z) {
    return new SWIGTYPE_p_std__arrayT_double_5_t(cStopPowJNI.AtomicData_get_shell_coeff(Z), true);
  }

  public AtomicData() {
    this(cStopPowJNI.new_AtomicData(), true);
  }

  public final static int n = cStopPowJNI.AtomicData_n_get();
}
