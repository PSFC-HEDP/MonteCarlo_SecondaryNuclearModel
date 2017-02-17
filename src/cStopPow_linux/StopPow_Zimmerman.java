/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 3.0.7
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

package cStopPow;

public class StopPow_Zimmerman extends StopPow_PartialIoniz {
  private transient long swigCPtr;

  protected StopPow_Zimmerman(long cPtr, boolean cMemoryOwn) {
    super(cStopPowJNI.StopPow_Zimmerman_SWIGUpcast(cPtr), cMemoryOwn);
    swigCPtr = cPtr;
  }

  protected static long getCPtr(StopPow_Zimmerman obj) {
    return (obj == null) ? 0 : obj.swigCPtr;
  }

  protected void finalize() {
    delete();
  }

  public synchronized void delete() {
    if (swigCPtr != 0) {
      if (swigCMemOwn) {
        swigCMemOwn = false;
        cStopPowJNI.delete_StopPow_Zimmerman(swigCPtr);
      }
      swigCPtr = 0;
    }
    super.delete();
  }

  public StopPow_Zimmerman(double mt, double Zt, DoubleVector mf, DoubleVector Zf, DoubleVector Tf, DoubleVector nf, DoubleVector Zbar, double Te) throws java.lang.IllegalArgumentException {
    this(cStopPowJNI.new_StopPow_Zimmerman__SWIG_0(mt, Zt, DoubleVector.getCPtr(mf), mf, DoubleVector.getCPtr(Zf), Zf, DoubleVector.getCPtr(Tf), Tf, DoubleVector.getCPtr(nf), nf, DoubleVector.getCPtr(Zbar), Zbar, Te), true);
  }

  public StopPow_Zimmerman(double mt, double Zt, SWIGTYPE_p_std__vectorT_std__arrayT_double_5_t_t field, double Te) throws java.lang.IllegalArgumentException {
    this(cStopPowJNI.new_StopPow_Zimmerman__SWIG_1(mt, Zt, SWIGTYPE_p_std__vectorT_std__arrayT_double_5_t_t.getCPtr(field), Te), true);
  }

  public double dEdx_MeV_um(double E) throws java.lang.IllegalArgumentException {
    return cStopPowJNI.StopPow_Zimmerman_dEdx_MeV_um(swigCPtr, this, E);
  }

  public double dEdx_MeV_mgcm2(double E) throws java.lang.IllegalArgumentException {
    return cStopPowJNI.StopPow_Zimmerman_dEdx_MeV_mgcm2(swigCPtr, this, E);
  }

  public double get_Emin() {
    return cStopPowJNI.StopPow_Zimmerman_get_Emin(swigCPtr, this);
  }

  public double get_Emax() {
    return cStopPowJNI.StopPow_Zimmerman_get_Emax(swigCPtr, this);
  }

  public double dEdx_free_electron(double E) {
    return cStopPowJNI.StopPow_Zimmerman_dEdx_free_electron(swigCPtr, this, E);
  }

  public double dEdx_bound_electron(double E) {
    return cStopPowJNI.StopPow_Zimmerman_dEdx_bound_electron(swigCPtr, this, E);
  }

  public double dEdx_ion(double E) {
    return cStopPowJNI.StopPow_Zimmerman_dEdx_ion(swigCPtr, this, E);
  }

  public void set_quantum(boolean set) {
    cStopPowJNI.StopPow_Zimmerman_set_quantum(swigCPtr, this, set);
  }

}
