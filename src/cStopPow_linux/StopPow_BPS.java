/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 3.0.7
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

package cStopPow;

public class StopPow_BPS extends StopPow_Plasma {
  private transient long swigCPtr;

  protected StopPow_BPS(long cPtr, boolean cMemoryOwn) {
    super(cStopPowJNI.StopPow_BPS_SWIGUpcast(cPtr), cMemoryOwn);
    swigCPtr = cPtr;
  }

  protected static long getCPtr(StopPow_BPS obj) {
    return (obj == null) ? 0 : obj.swigCPtr;
  }

  protected void finalize() {
    delete();
  }

  public synchronized void delete() {
    if (swigCPtr != 0) {
      if (swigCMemOwn) {
        swigCMemOwn = false;
        cStopPowJNI.delete_StopPow_BPS(swigCPtr);
      }
      swigCPtr = 0;
    }
    super.delete();
  }

  public StopPow_BPS(double mt, double Zt, DoubleVector mf, DoubleVector Zf, DoubleVector Tf, DoubleVector nf) throws java.lang.IllegalArgumentException {
    this(cStopPowJNI.new_StopPow_BPS__SWIG_0(mt, Zt, DoubleVector.getCPtr(mf), mf, DoubleVector.getCPtr(Zf), Zf, DoubleVector.getCPtr(Tf), Tf, DoubleVector.getCPtr(nf), nf), true);
  }

  public StopPow_BPS(double mt, double Zt, SWIGTYPE_p_std__vectorT_std__arrayT_double_4_t_t field) throws java.lang.IllegalArgumentException {
    this(cStopPowJNI.new_StopPow_BPS__SWIG_1(mt, Zt, SWIGTYPE_p_std__vectorT_std__arrayT_double_4_t_t.getCPtr(field)), true);
  }

  public StopPow_BPS(double mt, double Zt, DoubleVector mf, DoubleVector Zf, DoubleVector Tf, DoubleVector nf, double Te) throws java.lang.IllegalArgumentException {
    this(cStopPowJNI.new_StopPow_BPS__SWIG_2(mt, Zt, DoubleVector.getCPtr(mf), mf, DoubleVector.getCPtr(Zf), Zf, DoubleVector.getCPtr(Tf), Tf, DoubleVector.getCPtr(nf), nf, Te), true);
  }

  public StopPow_BPS(double mt, double Zt, SWIGTYPE_p_std__vectorT_std__arrayT_double_4_t_t field, double Te) throws java.lang.IllegalArgumentException {
    this(cStopPowJNI.new_StopPow_BPS__SWIG_3(mt, Zt, SWIGTYPE_p_std__vectorT_std__arrayT_double_4_t_t.getCPtr(field), Te), true);
  }

  public double dEdx_MeV_um(double E) throws java.lang.IllegalArgumentException {
    return cStopPowJNI.StopPow_BPS_dEdx_MeV_um(swigCPtr, this, E);
  }

  public double dEdx_MeV_mgcm2(double E) throws java.lang.IllegalArgumentException {
    return cStopPowJNI.StopPow_BPS_dEdx_MeV_mgcm2(swigCPtr, this, E);
  }

  public double dEdx_field(double E, int i) throws java.lang.IllegalArgumentException {
    return cStopPowJNI.StopPow_BPS_dEdx_field(swigCPtr, this, E, i);
  }

  public double dEdx_short(double E) {
    return cStopPowJNI.StopPow_BPS_dEdx_short__SWIG_0(swigCPtr, this, E);
  }

  public double dEdx_short(double E, int i) {
    return cStopPowJNI.StopPow_BPS_dEdx_short__SWIG_1(swigCPtr, this, E, i);
  }

  public double dEdx_long(double E) {
    return cStopPowJNI.StopPow_BPS_dEdx_long__SWIG_0(swigCPtr, this, E);
  }

  public double dEdx_long(double E, int i) {
    return cStopPowJNI.StopPow_BPS_dEdx_long__SWIG_1(swigCPtr, this, E, i);
  }

  public double dEdx_quantum(double E) {
    return cStopPowJNI.StopPow_BPS_dEdx_quantum__SWIG_0(swigCPtr, this, E);
  }

  public double dEdx_quantum(double E, int i) {
    return cStopPowJNI.StopPow_BPS_dEdx_quantum__SWIG_1(swigCPtr, this, E, i);
  }

  public double get_Emin() {
    return cStopPowJNI.StopPow_BPS_get_Emin(swigCPtr, this);
  }

  public double get_Emax() {
    return cStopPowJNI.StopPow_BPS_get_Emax(swigCPtr, this);
  }

  public double Fc_real(double u) {
    return cStopPowJNI.StopPow_BPS_Fc_real(swigCPtr, this, u);
  }

  public double Fc_imag(double u) {
    return cStopPowJNI.StopPow_BPS_Fc_imag(swigCPtr, this, u);
  }

}
