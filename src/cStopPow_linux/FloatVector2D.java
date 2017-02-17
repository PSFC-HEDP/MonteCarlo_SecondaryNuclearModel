/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 3.0.7
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

package cStopPow;

public class FloatVector2D {
  private transient long swigCPtr;
  protected transient boolean swigCMemOwn;

  protected FloatVector2D(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  protected static long getCPtr(FloatVector2D obj) {
    return (obj == null) ? 0 : obj.swigCPtr;
  }

  protected void finalize() {
    delete();
  }

  public synchronized void delete() {
    if (swigCPtr != 0) {
      if (swigCMemOwn) {
        swigCMemOwn = false;
        cStopPowJNI.delete_FloatVector2D(swigCPtr);
      }
      swigCPtr = 0;
    }
  }

  public FloatVector2D() {
    this(cStopPowJNI.new_FloatVector2D__SWIG_0(), true);
  }

  public FloatVector2D(long n) {
    this(cStopPowJNI.new_FloatVector2D__SWIG_1(n), true);
  }

  public long size() {
    return cStopPowJNI.FloatVector2D_size(swigCPtr, this);
  }

  public long capacity() {
    return cStopPowJNI.FloatVector2D_capacity(swigCPtr, this);
  }

  public void reserve(long n) {
    cStopPowJNI.FloatVector2D_reserve(swigCPtr, this, n);
  }

  public boolean isEmpty() {
    return cStopPowJNI.FloatVector2D_isEmpty(swigCPtr, this);
  }

  public void clear() {
    cStopPowJNI.FloatVector2D_clear(swigCPtr, this);
  }

  public void add(FloatVector x) {
    cStopPowJNI.FloatVector2D_add(swigCPtr, this, FloatVector.getCPtr(x), x);
  }

  public FloatVector get(int i) {
    return new FloatVector(cStopPowJNI.FloatVector2D_get(swigCPtr, this, i), false);
  }

  public void set(int i, FloatVector val) {
    cStopPowJNI.FloatVector2D_set(swigCPtr, this, i, FloatVector.getCPtr(val), val);
  }

}