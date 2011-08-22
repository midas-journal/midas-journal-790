/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImageBSplineCoefficients.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkImageBSplineCoefficients - convert image to b-spline knots
// .SECTION Description
// vtkImageBSplineCoefficients prepares an image for b-spline
// interpolation by converting the image values into b-spline
// knot coefficients.  It is a necessary pre-filtering step
// before applying b-spline interpolation with vtkImageReslice.
//
// This class is based on code provided by Philippe Thevenaz of
// EPFL, Lausanne, Switzerland.  Please acknowledge his contribution
// by citing the following paper:
// [1] P. Thevenaz, T. Blu, M. Unser, "Interpolation Revisited,"
//     IEEE Transactions on Medical Imaging 19(7):739-758, 2000.
//
// The clamped boundary condition (which is the default) is taken
// from code presented in the following paper:
// [2] D. Ruijters, P. Thevenaz,
//     "GPU Prefilter for Accurate Cubic B-spline Interpolation,"
//     The Computer Journal, doi: 10.1093/comjnl/bxq086, 2010.
//
// .SECTION Thanks
// This class was written by David Gobbi at the Seaman Family MR Research
// Centre, Foothills Medical Centre, Calgary, Alberta.


#ifndef __vtkImageBSplineCoefficients_h
#define __vtkImageBSplineCoefficients_h


#include "vtkThreadedImageAlgorithm.h"

#define VTK_BSPLINE_CLAMP 0
#define VTK_BSPLINE_MIRROR  1
#define VTK_BSPLINE_REPEAT  2

class VTK_EXPORT vtkImageBSplineCoefficients :
  public vtkThreadedImageAlgorithm
{
public:
  static vtkImageBSplineCoefficients *New();
  vtkTypeMacro(vtkImageBSplineCoefficients,vtkThreadedImageAlgorithm);

  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Set the degree of the spline polynomial.  The default value is 3,
  // and the maximum is 9.
  vtkSetClampMacro(SplineDegree, int, 0, 9);
  vtkGetMacro(SplineDegree, int);

  // Description:
  // Set the border mode.  The filter that is used to create the
  // coefficients must repeat the image somehow to make a theoritically
  // infinite input.  The default is to clamp values that are off the
  // edge of the image, to the value at the closest point on the edge.
  // The other ways of virtually extending the image are to produce
  // mirrored copies, which results in optimal smoothness at the boundary,
  // or to repeat the image, which results in a cyclic or periodic spline.
  vtkSetClampMacro(BorderMode, int, VTK_BSPLINE_CLAMP, VTK_BSPLINE_REPEAT);
  void SetBorderModeToClamp() { this->SetBorderMode(VTK_BSPLINE_CLAMP); }
  void SetBorderModeToMirror() { this->SetBorderMode(VTK_BSPLINE_MIRROR); }
  void SetBorderModeToRepeat() { this->SetBorderMode(VTK_BSPLINE_REPEAT); }
  vtkGetMacro(BorderMode, int);
  const char *GetBorderModeAsString();

  // Description:
  // Set the scalar type of the output.  Default is float.
  // Floating-point output is used to avoid overflow, since the
  // range of the output values is larger than the input values.
  vtkSetClampMacro(OutputScalarType, int, VTK_FLOAT, VTK_DOUBLE);
  vtkGetMacro(OutputScalarType, int);
  void SetOutputScalarTypeToFloat() {
    this->SetOutputScalarType(VTK_FLOAT); }
  void SetOutputScalarTypeToDouble() {
    this->SetOutputScalarType(VTK_DOUBLE); }
  const char *GetOutputScalarTypeAsString();

  // Description:
  // Bypass the filter, do not do any processing.  If this is on,
  // then the output data will reference the input data directly,
  // and the output type will be the same as the input type.  This
  // is useful a downstream filter sometimes uses b-spline interpolation
  // and sometimes uses other forms of interpolation.
  vtkSetMacro(Bypass, int);
  vtkBooleanMacro(Bypass, int);
  vtkGetMacro(Bypass, int);

  // Description:
  // Check a point against the image bounds.  Return 0 if out of bounds,
  // and 1 if inside bounds.  Calling Evaluate on a point outside the
  // bounds will not generate an error, but the value returned will
  // depend on the BorderMode.
  int CheckBounds(const double point[3]);

  // Description:
  // Interpolate a value from the image.  You must call Update() before
  // calling this method for the first time.  The first signature can
  // return multiple components, while the second signature is for use
  // on single-component images.
  void Evaluate(const double point[3], double *value);
  double Evaluate(double x, double y, double z);
  double Evaluate(const double point[3]) {
    return this->Evaluate(point[0], point[1], point[2]); }

  // Description:
  // Internal method.  Override SplitExtent so that the full extent is
  // available in the direction currently being processed.
  int SplitExtent(int splitExt[6], int startExt[6], int num, int total);

  // Description:
  // Internal method.  Get the poles for spline of given degree.
  // Returns zero if an illegal degree is given (allowed range 2 to 9).
  // The parameter numPoles will be set to a value between 1 and 4.
  static int GetPoleValues(double poles[4], long &numPoles, long degree);

  // Description:
  // Internal method.  Compute the coefficients for one row of data.
  static void ConvertToInterpolationCoefficients(
    double data[], long size, long border, double poles[4], long numPoles,
    double tol);

  // Description:
  // Internal method.  Get interpolation weights for offset w, where
  // w is between 0 and 1.  You must provide the degree of the spline.
  static int GetInterpolationWeights(
    double weights[10], double w, long degree);

protected:
  vtkImageBSplineCoefficients();
  ~vtkImageBSplineCoefficients();

  virtual void AllocateOutputData(vtkImageData *out, int *uExtent);
  virtual vtkImageData *AllocateOutputData(vtkDataObject *out);

  virtual int RequestData(
    vtkInformation*, vtkInformationVector**, vtkInformationVector*);
  virtual int RequestInformation(
    vtkInformation*, vtkInformationVector**, vtkInformationVector*);
  virtual int RequestUpdateExtent(
    vtkInformation*, vtkInformationVector**, vtkInformationVector*);

  virtual void ThreadedExecute(vtkImageData *inData, vtkImageData *outData,
                               int outExt[6], int threadId);

  static double InitialCausalCoefficient(
    double data[], long size, long border, double pole, double tol);

  static double InitialAntiCausalCoefficient(
    double data[], long size, long border, double pole, double tol);

//BTX
  template<class T>
  static int InterpolatedValue(
    const T *coeffs, T *value,
    long width, long height, long slices, long depth,
    double x, double y, double z, long degree, long border);
//ETX

  int SplineDegree;
  int BorderMode;
  int OutputScalarType;
  int Bypass;
  int DataWasPassed;
  int Iteration;

private:
  vtkImageBSplineCoefficients(const vtkImageBSplineCoefficients&);  // Not implemented.
  void operator=(const vtkImageBSplineCoefficients&);  // Not implemented.
};

#endif
