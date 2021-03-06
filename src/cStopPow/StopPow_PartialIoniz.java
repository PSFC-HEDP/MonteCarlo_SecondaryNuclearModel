/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 3.0.10
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

package cStopPow;

public class StopPow_PartialIoniz extends StopPow {
  private transient long swigCPtr;

  protected StopPow_PartialIoniz(long cPtr, boolean cMemoryOwn) {
    super(cStopPowJNI.StopPow_PartialIoniz_SWIGUpcast(cPtr), cMemoryOwn);
    swigCPtr = cPtr;
  }

  protected static long getCPtr(StopPow_PartialIoniz obj) {
    return (obj == null) ? 0 : obj.swigCPtr;
  }

  protected void finalize() {
    delete();
  }

  public synchronized void delete() {
    if (swigCPtr != 0) {
      if (swigCMemOwn) {
        swigCMemOwn = false;
        cStopPowJNI.delete_StopPow_PartialIoniz(swigCPtr);
      }
      swigCPtr = 0;
    }
    super.delete();
  }

  public void set_particle(double mt, double Zt) {
    cStopPowJNI.StopPow_PartialIoniz_set_particle(swigCPtr, this, mt, Zt);
  }

  public void set_field(DoubleVector mf, DoubleVector Zf, DoubleVector Tf, DoubleVector nf, DoubleVector Zbar, double Te) {
    cStopPowJNI.StopPow_PartialIoniz_set_field__SWIG_0(swigCPtr, this, DoubleVector.getCPtr(mf), mf, DoubleVector.getCPtr(Zf), Zf, DoubleVector.getCPtr(Tf), Tf, DoubleVector.getCPtr(nf), nf, DoubleVector.getCPtr(Zbar), Zbar, Te);
  }

  public void set_field(SWIGTYPE_p_std__vectorT_std__arrayT_double_5_t_t field, double Te) throws java.lang.IllegalArgumentException {
    cStopPowJNI.StopPow_PartialIoniz_set_field__SWIG_1(swigCPtr, this, SWIGTYPE_p_std__vectorT_std__arrayT_double_5_t_t.getCPtr(field), Te);
  }

}
