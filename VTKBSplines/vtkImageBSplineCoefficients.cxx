/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImageBSplineCoefficients.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkImageBSplineCoefficients.h"

#include "vtkMath.h"
#include "vtkIntArray.h"
#include "vtkImageData.h"
#include "vtkPointData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkTemplateAliasMacro.h"

#include <stddef.h>

vtkStandardNewMacro(vtkImageBSplineCoefficients);

//----------------------------------------------------------------------------
vtkImageBSplineCoefficients::vtkImageBSplineCoefficients()
{
  this->SplineDegree = 3;
  this->BorderMode = VTK_BSPLINE_CLAMP;
  this->OutputScalarType = VTK_FLOAT;
  this->Bypass = 0;
  this->DataWasPassed = 0;
  this->Iteration = 0;
}

//----------------------------------------------------------------------------
vtkImageBSplineCoefficients::~vtkImageBSplineCoefficients()
{
}

//----------------------------------------------------------------------------
void vtkImageBSplineCoefficients::AllocateOutputData(
  vtkImageData *vtkNotUsed(output), int *vtkNotUsed(uExtent))
{
  // turn into a no-op, we allocate our output manually
}

//----------------------------------------------------------------------------
vtkImageData *vtkImageBSplineCoefficients::AllocateOutputData(
  vtkDataObject *output)
{
  // turn into a no-op, we allocate our output manually
  vtkImageData *out = vtkImageData::SafeDownCast(output);
  return out;
}

//----------------------------------------------------------------------------
int vtkImageBSplineCoefficients::RequestData(
  vtkInformation* request, vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
  vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation* outInfo = outputVector->GetInformationObject(0);

  vtkImageData *inData = vtkImageData::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkImageData *outData = vtkImageData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  outData->GetFieldData()->RemoveArray("BSpline");

  if (this->Bypass)
    {
    // directly pass the scalars to the output
    outData->SetExtent(inData->GetExtent());
    outData->GetPointData()->PassData(inData->GetPointData());
    this->DataWasPassed = 1;
    return 1;
    }
  else if (this->DataWasPassed)
    {
    // force reallocation of the scalars
    outData->GetPointData()->SetScalars(NULL);
    this->DataWasPassed = 0;
    }

  // Allocate the output data
  outData->SetExtent(outData->GetWholeExtent());
  outData->AllocateScalars();

  if (outData->GetScalarType() != VTK_FLOAT &&
      outData->GetScalarType() != VTK_DOUBLE)
    {
    vtkErrorMacro(<< "Execute: output data must be be type float or double.");
    return 0;
    }

  // copy the data to the output
  int ie[6], oe[6];
  inData->GetExtent(ie);
  outData->GetExtent(oe);
  if (ie[0] == oe[0] && ie[1] == oe[1] && ie[2] == oe[2] &&
      ie[3] == oe[3] && ie[4] == oe[4] && ie[5] == oe[5])
    {
    outData->GetPointData()->GetScalars()->DeepCopy(
      inData->GetPointData()->GetScalars());
    }
  else
    {
    vtkErrorMacro(<< "Execute: input and output extents do not match: "
                  << "(" << ie[0] << "," << ie[1] << "," << ie[2] << ","
                  << ie[3] << "," << ie[4] << "," << ie[5] << ") vs. "
                  << "(" << oe[0] << "," << oe[1] << "," << oe[2] << ","
                  << oe[3] << "," << oe[4] << "," << oe[5] << ")");
    return 0;
    }

  // Add a field data array with bspline information
  vtkIntArray *info = vtkIntArray::New();
  info->SetName("BSpline");
  info->SetNumberOfValues(2);
  info->SetValue(0, this->SplineDegree);
  info->SetValue(1, this->BorderMode);
  outData->GetFieldData()->AddArray(info);
  info->Delete();

  // if spline degree is < 2, no operation is required
  if (this->SplineDegree < 2)
    {
    return 1;
    }

  // execute over the three directions
  for (int i = 0; i < 3; i++)
    {
    this->Iteration = i;

    if (ie[2*i+1] > ie[2*i])
      {
      if (!this->vtkThreadedImageAlgorithm::RequestData(
            request, &outputVector, outputVector))
        {
        return 0;
        }
      }
    }

  return 1;
}

//----------------------------------------------------------------------------
int vtkImageBSplineCoefficients::RequestInformation(
  vtkInformation* vtkNotUsed(request),
  vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  int numComponents = 1;
  int scalarType = VTK_FLOAT;

  vtkInformation *inScalarInfo =
    vtkDataObject::GetActiveFieldInformation(inInfo,
      vtkDataObject::FIELD_ASSOCIATION_POINTS,
      vtkDataSetAttributes::SCALARS);

  if (inScalarInfo)
    {
    if (inScalarInfo->Has(vtkDataObject::FIELD_NUMBER_OF_COMPONENTS()))
      {
      numComponents =
        inScalarInfo->Get(vtkDataObject::FIELD_NUMBER_OF_COMPONENTS());
      }
    scalarType = inScalarInfo->Get(vtkDataObject::FIELD_ARRAY_TYPE());
    }

  if (this->Bypass)
    {
    vtkDataObject::SetPointDataActiveScalarInfo(
      outInfo, scalarType, numComponents);
    }
  else if (this->OutputScalarType == VTK_DOUBLE)
    {
    vtkDataObject::SetPointDataActiveScalarInfo(
      outInfo, VTK_DOUBLE, numComponents);
    }
  else
    {
    vtkDataObject::SetPointDataActiveScalarInfo(
      outInfo, VTK_FLOAT, numComponents);
    }

  return 1;
}

//----------------------------------------------------------------------------
int vtkImageBSplineCoefficients::RequestUpdateExtent(
  vtkInformation* vtkNotUsed(request),
  vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  int extent[6];

  if (this->Bypass)
    {
    // in bypass mode, just pass the update extent
    outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), extent);
    }
  else
    {
    // the whole input extent is required every time
    inInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), extent);
    }

  inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), extent, 6);

  return 1;
}

//----------------------------------------------------------------------------
// For streaming and threads.  Splits output update extent into num pieces.
// This method needs to be called num times.  Results must not overlap for
// consistent starting extent.  This particular filter requires splitting
// along different directions depending on the iteration.
int vtkImageBSplineCoefficients::SplitExtent(
  int splitExt[6], int fullExt[6], int num, int total)
{
  static int splitAxisPriority[3][2] = {
    { 2, 1 }, { 2, 0 }, { 1, 0 } };

  // start with same extent
  for (int i = 0; i < 6; i++)
    {
    splitExt[i] = fullExt[i];
    }

  int axis = this->Iteration;
  int splitAxis = splitAxisPriority[axis][0];
  int minIdx = fullExt[splitAxis*2];
  int maxIdx = fullExt[splitAxis*2 + 1];
  int size = maxIdx - minIdx + 1;
  if (size == 1)
    {
    splitAxis = splitAxisPriority[axis][1];
    minIdx = fullExt[splitAxis*2];
    maxIdx = fullExt[splitAxis*2 + 1];
    size = maxIdx - minIdx + 1;
    }

  // determine the actual number of pieces that will be generated
  if (size < total)
    {
    total = size;
    }

  // make sure that num isn't greater than the number of possible splits
  if (num < total)
    {
    splitExt[splitAxis*2] = minIdx + size*num/total;
    splitExt[splitAxis*2 + 1] = minIdx + size*(num + 1)/total - 1;
    }

  return total;
}

//----------------------------------------------------------------------------
template <class T>
void vtkImageBSplineCoefficientsExecute(
  vtkImageBSplineCoefficients* self,
  vtkImageData* inData, vtkImageData* outData, T *inPtr, T *outPtr,
  int extent[6], int axis, int threadId)
{
  // change the order so the inner loop is the chosen axis
  static int permute[3][3] = {
    { 0, 1, 2 }, { 1, 0, 2 }, { 2, 0, 1 } };

  int borderMode = self->GetBorderMode();

  int inExtent[6];
  inData->GetExtent(inExtent);
  int inMin0 = inExtent[2*permute[axis][0]];
  int inMax0 = inExtent[2*permute[axis][0] + 1];

  int outMin0 = extent[2*permute[axis][0]];
  int outMax0 = extent[2*permute[axis][0] + 1];
  int outMin1 = extent[2*permute[axis][1]];
  int outMax1 = extent[2*permute[axis][1] + 1];
  int outMin2 = extent[2*permute[axis][2]];
  int outMax2 = extent[2*permute[axis][2] + 1];

  vtkIdType inInc[6];
  inData->GetIncrements(inInc);
  vtkIdType inInc0 = inInc[permute[axis][0]];
  vtkIdType inInc1 = inInc[permute[axis][1]];
  vtkIdType inInc2 = inInc[permute[axis][2]];

  vtkIdType outInc[6];
  outData->GetIncrements(outInc);
  vtkIdType outInc0 = outInc[permute[axis][0]];
  vtkIdType outInc1 = outInc[permute[axis][1]];
  vtkIdType outInc2 = outInc[permute[axis][2]];

  int numscalars = outData->GetNumberOfScalarComponents();

  // for progress reporting
  unsigned long count = 0;
  unsigned long target = static_cast<unsigned long>(
    0.02*(outMax2-outMin2+1)*(outMax1-outMin1+1));
  target++;

  // Get the poles for the spline
  double poles[4];
  long numPoles;
  self->GetPoleValues(poles, numPoles, self->GetSplineDegree());

  // allocate workspace for one row
  double* image = new double[inMax0 - inMin0 + 1];

  // loop over all the extra axes
  T *inPtr2 = inPtr - (outMin0 - inMin0)*inInc0;
  T *outPtr2 = outPtr;
  for (int idx2 = outMin2; idx2 <= outMax2; idx2++)
    {
    T *inPtr1 = inPtr2;
    T *outPtr1 = outPtr2;
    for (int idx1 = outMin1; !self->AbortExecute && idx1 <= outMax1; idx1++)
      {
      if (threadId == 0 && count % target == 0)
        {
        self->UpdateProgress((axis + count/(50.0*target))/3.0);
        }
      count++;

      // loop over components
      for (int idxC = 0; idxC < numscalars; idxC++)
        {
        T *inPtr0 = inPtr1 + idxC;
        T *outPtr0 = outPtr1 + idxC;

        double *imagePtr = image;
        for (int jdx0 = inMin0; jdx0 <= inMax0; jdx0++)
          {
          *imagePtr++ = static_cast<double>(*inPtr0);
          inPtr0 += inInc0;
          }

        // Call the code that generates the b-spline knots,
        self->ConvertToInterpolationCoefficients(
          image, inMax0 - inMin0 + 1, borderMode, poles, numPoles,
          VTK_DBL_EPSILON);

        // Copy to output
        imagePtr = image + (outMin0 - inMin0);
        for (int idx0 = outMin0; idx0 <= outMax0; idx0++)
          {
          *outPtr0 = *imagePtr++;
          outPtr0 += outInc0;
          }
        }
      inPtr1 += inInc1;
      outPtr1 += outInc1;
      }
    inPtr2 += inInc2;
    outPtr2 += outInc2;
    }

  delete [] image;
}

//----------------------------------------------------------------------------
// This is called three times (once per dimension)
void vtkImageBSplineCoefficients::ThreadedExecute(
  vtkImageData *inData, vtkImageData *outData, int outExt[6], int threadId)
{
  void *inPtr = inData->GetScalarPointerForExtent(outExt);
  void *outPtr = outData->GetScalarPointerForExtent(outExt);

  if (outData->GetScalarType() == VTK_FLOAT)
    {
    vtkImageBSplineCoefficientsExecute(
      this, inData, outData,
      static_cast<float*>(inPtr), static_cast<float*>(outPtr),
      outExt, this->Iteration, threadId);
    }
  else if (outData->GetScalarType() == VTK_DOUBLE)
    {
    vtkImageBSplineCoefficientsExecute(
      this, inData, outData,
      static_cast<double*>(inPtr), static_cast<double*>(outPtr),
      outExt, this->Iteration, threadId);
    }
}

//----------------------------------------------------------------------------
void vtkImageBSplineCoefficients::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  os << "SplineDegree: " << this->SplineDegree << "\n";
  os << "BorderMode: " << this->GetBorderModeAsString() << "\n";
  os << "OutputScalarType: " << this->GetOutputScalarTypeAsString() << "\n";
  os << "Bypass: " << (this->Bypass ? "On\n" : "Off\n" );
}

//----------------------------------------------------------------------------
const char *vtkImageBSplineCoefficients::GetBorderModeAsString()
{
  switch (this->BorderMode)
    {
    case VTK_BSPLINE_CLAMP:
      return "Clamp";
    case VTK_BSPLINE_MIRROR:
      return "Mirror";
    case VTK_BSPLINE_REPEAT:
      return "Repeat";
    default:
      break;
    }

  return "Unknown";
}

//----------------------------------------------------------------------------
const char *vtkImageBSplineCoefficients::GetOutputScalarTypeAsString()
{
  return vtkImageScalarTypeNameMacro(this->OutputScalarType);
}

//----------------------------------------------------------------------------
int vtkImageBSplineCoefficients::CheckBounds(const double point[3])
{
  double *bounds = this->GetOutput()->GetBounds();
  for (int i = 0; i < 3; i++)
    {
    double a = bounds[0];
    double b = bounds[1];
    if ((b - a) > 1e-16 && (point[i] < a || point[i] > b))
      {
      return 0;
      }
    bounds++;
    }
  return 1;
}

//----------------------------------------------------------------------------
void vtkImageBSplineCoefficients::Evaluate(const double p[3], double *val)
{
  vtkImageData *output = this->GetOutput();
  int extent[6];
  output->GetExtent(extent);
  int width = extent[1] - extent[0] + 1;
  int height = extent[3] - extent[2] + 1;
  int slices = extent[5] - extent[4] + 1;

  if (width < 1 || height < 1 || slices < 1)
    {
    vtkErrorMacro("Called Evaluate() on empty extent");
    return;
    }

  int numscalars = output->GetNumberOfScalarComponents();
  int scalarType = output->GetScalarType();

  if (scalarType == VTK_FLOAT)
    {
    float *coeffs = static_cast<float *>(output->GetScalarPointer());
    float value4[4];
    float *value = value4;
    if (numscalars > 4)
      {
      value = new float[numscalars];
      }

    this->InterpolatedValue(
      coeffs, value, width, height, slices, numscalars, p[0], p[1], p[2],
      this->SplineDegree, this->BorderMode);

    for (int i = 0; i < numscalars; i++)
      {
      val[i] = value[i];
      }

    if (value != value4)
      {
      delete [] value;
      }
    }
  else if (scalarType == VTK_DOUBLE)
    {
    double *coeffs = static_cast<double *>(output->GetScalarPointer());
    this->InterpolatedValue(
      coeffs, val, width, height, slices, numscalars, p[0], p[1], p[2],
      this->SplineDegree, this->BorderMode);
    }
  else
    {
    vtkErrorMacro("Called Evaluate(), but data is not float or double.");
    }
}

//----------------------------------------------------------------------------
double vtkImageBSplineCoefficients::Evaluate(double x, double y, double z)
{
  vtkImageData *output = this->GetOutput();
  int extent[6];
  output->GetExtent(extent);
  int width = extent[1] - extent[0] + 1;
  int height = extent[3] - extent[2] + 1;
  int slices = extent[5] - extent[4] + 1;

  if (width < 1 || height < 1 || slices < 1)
    {
    vtkErrorMacro("Called Evaluate() on empty extent");
    return 0.0;
    }

  int numscalars = output->GetNumberOfScalarComponents();
  int scalarType = output->GetScalarType();

  if (scalarType == VTK_FLOAT)
    {
    float *coeffs = static_cast<float *>(output->GetScalarPointer());
    float value4[4];
    float *value = value4;
    if (numscalars > 4)
      {
      value = new float[numscalars];
      }

    this->InterpolatedValue(
      coeffs, value, width, height, slices, numscalars, x, y, z,
      this->SplineDegree, this->BorderMode);

    if (value != value4)
      {
      value4[0] = value[0];
      delete [] value;
      value = value4;
      }

    return value[0];
    }
  else if (scalarType == VTK_DOUBLE)
    {
    double *coeffs = static_cast<double *>(output->GetScalarPointer());
    double value4[4];
    double *value = value4;
    if (numscalars > 4)
      {
      value = new double[numscalars];
      }

    this->InterpolatedValue(
      coeffs, value, width, height, slices, numscalars, x, y, z,
      this->SplineDegree, this->BorderMode);

    if (value != value4)
      {
      value4[0] = value[0];
      delete [] value;
      value = value4;
      }

    return value[0];
    }
  else
    {
    vtkErrorMacro("Called Evaluate(), but data is not float or double.");
    }

  return 0;
}


/*****************************************************************************
 *    Date: January 5, 2009
 *----------------------------------------------------------------------------
 *    This C program is based on the following paper:
 *        P. Thevenaz, T. Blu, M. Unser, "Interpolation Revisited,"
 *        IEEE Transactions on Medical Imaging,
 *        vol. 19, no. 7, pp. 739-758, July 2000.
 *----------------------------------------------------------------------------
 *    Philippe Thevenaz
 *    EPFL/STI/IMT/LIB/BM.4.137
 *    Station 17
 *    CH-1015 Lausanne VD
 *----------------------------------------------------------------------------
 *    phone (CET):    +41(21)693.51.61
 *    fax:            +41(21)693.37.01
 *    RFC-822:        philippe.thevenaz@epfl.ch
 *    X-400:            /C=ch/A=400net/P=switch/O=epfl/S=thevenaz/G=philippe/
 *    URL:            http://bigwww.epfl.ch/
 ****************************************************************************/

/*--------------------------------------------------------------------------*/
void vtkImageBSplineCoefficients::
ConvertToInterpolationCoefficients
    (
        double  c[],        /* input samples --> output coefficients */
        long    DataLength, /* number of samples or coefficients */
        long    Border,     /* border mode */
        double  z[],        /* poles */
        long    NbPoles,    /* number of poles */
        double  Tolerance   /* admissible relative error */
    )

{ /* begin ConvertToInterpolationCoefficients */

    double  Lambda = 1.0;
    long    n, k;

    /* special case required by mirror boundaries */
    if (DataLength == 1L) {
        return;
    }
    /* compute the overall gain */
    for (k = 0L; k < NbPoles; k++) {
        Lambda = Lambda * (1.0 - z[k]) * (1.0 - 1.0 / z[k]);
    }
    /* apply the gain */
    for (n = 0L; n < DataLength; n++) {
        c[n] *= Lambda;
    }
    /* loop over all poles */
    for (k = 0L; k < NbPoles; k++) {
        /* causal initialization */
        c[0] = InitialCausalCoefficient(c, DataLength, Border, z[k], Tolerance);
        /* causal recursion */
        for (n = 1L; n < DataLength; n++) {
            c[n] += z[k] * c[n - 1L];
        }
        /* anticausal initialization */
        c[DataLength - 1L] = InitialAntiCausalCoefficient(c, DataLength, Border, z[k], Tolerance);
        /* anticausal recursion */
        for (n = DataLength - 2L; 0 <= n; n--) {
            c[n] = z[k] * (c[n + 1L] - c[n]);
        }
    }
} /* end ConvertToInterpolationCoefficients */

/*--------------------------------------------------------------------------*/
double vtkImageBSplineCoefficients::
InitialCausalCoefficient
    (
        double  c[],        /* coefficients */
        long    DataLength, /* number of coefficients */
        long    Border,     /* border mode */
        double  z,          /* actual pole */
        double  Tolerance   /* admissible relative error */
    )

{ /* begin InitialCausalCoefficient */

    double  Sum, zn, z2n, iz;
    long    n, Horizon;

    switch (Border) {
        case VTK_BSPLINE_CLAMP:
            /* this initialization corresponds to repeating edge pixels */
            Horizon = DataLength;
            if (Tolerance > 0.0 && DataLength > 16) {
                Horizon = (long)(ceil(log(Tolerance) / log(fabs(z))));
            }
            if (Horizon < DataLength) {
                /* accelerated loop */
                zn = z;
                Sum = c[0];
                for (n = 0; n < Horizon; n++) {
                    Sum += zn * c[n];
                    zn *= z;
                }
                return(Sum);
            }
            else {
                /* full loop */
                zn = z;
                iz = 1.0 / z;
                z2n = pow(z, (double)DataLength);
                Sum = z * c[0] + z2n * z2n * c[0];
                z2n *= z2n * iz;
                for (n = 1L; n <= DataLength - 1L; n++) {
                    zn *= z;
                    Sum += (zn + z2n) * c[n];
                    z2n *= iz;
                }
                return(c[0] + Sum / (1.0 - zn * zn));
            }
            break;

        case VTK_BSPLINE_MIRROR:
            /* this initialization corresponds to mirror boundaries */
            Horizon = DataLength;
            if (Tolerance > 0.0 && DataLength > 16) {
                Horizon = (long)(ceil(log(Tolerance) / log(fabs(z))));
            }
            if (Horizon < DataLength) {
                /* accelerated loop */
                zn = z;
                Sum = c[0];
                for (n = 1L; n < Horizon; n++) {
                    Sum += zn * c[n];
                    zn *= z;
                }
                return(Sum);
            }
            else {
                /* full loop */
                zn = z;
                iz = 1.0 / z;
                z2n = pow(z, (double)(DataLength - 1L));
                Sum = c[0] + z2n * c[DataLength - 1L];
                z2n *= z2n * iz;
                for (n = 1L; n <= DataLength - 2L; n++) {
                    Sum += (zn + z2n) * c[n];
                    zn *= z;
                    z2n *= iz;
                }
                return(Sum / (1.0 - zn * zn));
            }
            break;

        case VTK_BSPLINE_REPEAT:
            /* this initialization corresponds to periodic boundaries */
            Horizon = DataLength;
            if (Tolerance > 0.0 && DataLength > 16) {
                Horizon = (long)(ceil(log(Tolerance) / log(fabs(z))));
            }
            if (Horizon < DataLength) {
                /* accelerated loop */
                zn = z;
                Sum = c[0];
                for (n = 1L; n < Horizon; n++) {
                     Sum += zn * c[DataLength-n];
                     zn *= z;
                }
                return(Sum);
            }
            else {
                /* full loop */
                zn = z;
                Sum = c[0];
                for (n = 1L; n < DataLength; n++) {
                     Sum += zn * c[DataLength-n];
                     zn *= z;
                }
                return(Sum / (1.0 - zn));
            }
            break;
    }

    return(0.0);
} /* end InitialCausalCoefficient */

/*--------------------------------------------------------------------------*/
double vtkImageBSplineCoefficients::
InitialAntiCausalCoefficient
    (
        double  c[],        /* coefficients */
        long    DataLength, /* number of samples or coefficients */
        long    Border,     /* border mode */
        double  z,          /* actual pole */
        double  Tolerance   /* admissible relative error */
    )

{ /* begin InitialAntiCausalCoefficient */
    double Sum;
    double zn;
    long    n, Horizon;

    switch (Border) {
        case VTK_BSPLINE_CLAMP:
            /* this initialization corresponds to repeating edge pixels */
            return((z / (z - 1.0)) * c[DataLength - 1L]);
        case VTK_BSPLINE_MIRROR:
            /* this initialization corresponds to mirror boundaries */
            return((z / (z * z - 1.0)) * (z * c[DataLength - 2L] + c[DataLength - 1L]));
        case VTK_BSPLINE_REPEAT:
            /* this initialization corresponds to periodic boundaries */
            Horizon = DataLength;
            if (Tolerance > 0.0 && DataLength > 16) {
                Horizon = (long)(ceil(log(Tolerance) / log(fabs(z))));
            }
            if (Horizon < DataLength) {
                /* accelerated loop */
                zn = z;
                Sum = c[0];
                for (n = 1L; n < Horizon; n++) {
                    Sum += zn * c[n];
                    zn *= z;
                }
                return(-Sum * z * z - z * c[DataLength - 1L]);
            }
            else {
                /* full loop */
                zn = z;
                Sum = c[0];
                for (n = 1L; n < DataLength; n++) {
                     Sum += zn * c[n];
                     zn *= z;
                }
                return(Sum * z * z / (zn - 1.0) - z * c[DataLength - 1L]);
            }
    }
    return(0.0);
} /* end InitialAntiCausalCoefficient */

/*--------------------------------------------------------------------------*/
int vtkImageBSplineCoefficients::
GetPoleValues
    (
        double  Pole[4],     /* poles for b-spline */
        long   &NbPoles,     /* number of poles (will be four, at maximum) */
        long    SplineDegree  /* degree of the spline model */
    )

{ /* begin GetPoleValues */

    /* recover the poles from a lookup table */
    switch (SplineDegree) {
        case 0L:
            NbPoles = 0L;
            break;
        case 1L:
            NbPoles = 0L;
            break;
        case 2L:
            NbPoles = 1L;
            Pole[0] = sqrt(8.0) - 3.0;
            break;
        case 3L:
            NbPoles = 1L;
            Pole[0] = sqrt(3.0) - 2.0;
            break;
        case 4L:
            NbPoles = 2L;
            Pole[0] = sqrt(664.0 - sqrt(438976.0)) + sqrt(304.0) - 19.0;
            Pole[1] = sqrt(664.0 + sqrt(438976.0)) - sqrt(304.0) - 19.0;
            break;
        case 5L:
            NbPoles = 2L;
            Pole[0] = sqrt(135.0 / 2.0 - sqrt(17745.0 / 4.0)) + sqrt(105.0 / 4.0)
                - 13.0 / 2.0;
            Pole[1] = sqrt(135.0 / 2.0 + sqrt(17745.0 / 4.0)) - sqrt(105.0 / 4.0)
                - 13.0 / 2.0;
            break;
        case 6L:
            NbPoles = 3L;
            Pole[0] = -0.48829458930304475513011803888378906211227916123938;
            Pole[1] = -0.081679271076237512597937765737059080653379610398148;
            Pole[2] = -0.0014141518083258177510872439765585925278641690553467;
            break;
        case 7L:
            NbPoles = 3L;
            Pole[0] = -0.53528043079643816554240378168164607183392315234269;
            Pole[1] = -0.12255461519232669051527226435935734360548654942730;
            Pole[2] = -0.0091486948096082769285930216516478534156925639545994;
            break;
        case 8L:
            NbPoles = 4L;
            Pole[0] = -0.57468690924876543053013930412874542429066157804125;
            Pole[1] = -0.16303526929728093524055189686073705223476814550830;
            Pole[2] = -0.023632294694844850023403919296361320612665920854629;
            Pole[3] = -0.00015382131064169091173935253018402160762964054070043;
            break;
        case 9L:
            NbPoles = 4L;
            Pole[0] = -0.60799738916862577900772082395428976943963471853991;
            Pole[1] = -0.20175052019315323879606468505597043468089886575747;
            Pole[2] = -0.043222608540481752133321142979429688265852380231497;
            Pole[3] = -0.0021213069031808184203048965578486234220548560988624;
            break;
        default:
            NbPoles = 0L;
            return(1);
            break;
    }

    return(0);
} /* end GetPoleValues */

/*--------------------------------------------------------------------------*/
int vtkImageBSplineCoefficients::
GetInterpolationWeights
    (
        double xWeight[10],  /* weights, size is SplineDegree + 1 */
        double w,            /* offset, value between 0 and 1 */
        long SplineDegree    /* degree of spline, between 2 and 9 */
    )
{ /* begin GetInterpolationWeights */
    double    w2, w4, t, t0, t1;

    /* compute the interpolation weights */
    switch (SplineDegree) {
        case 0L:
            xWeight[0] = 1.0;
            break;
        case 1L:
            xWeight[0] = 1.0 - w;
            xWeight[1] = w;
            break;
        case 2L:
            xWeight[1] = 3.0 / 4.0 - w * w;
            xWeight[2] = (1.0 / 2.0) * (w - xWeight[1] + 1.0);
            xWeight[0] = 1.0 - xWeight[1] - xWeight[2];
            break;
        case 3L:
            xWeight[3] = (1.0 / 6.0) * w * w * w;
            xWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - xWeight[3];
            xWeight[2] = w + xWeight[0] - 2.0 * xWeight[3];
            xWeight[1] = 1.0 - xWeight[0] - xWeight[2] - xWeight[3];
            break;
        case 4L:
            w2 = w * w;
            t = (1.0 / 6.0) * w2;
            xWeight[0] = 1.0 / 2.0 - w;
            xWeight[0] *= xWeight[0];
            xWeight[0] *= (1.0 / 24.0) * xWeight[0];
            t0 = w * (t - 11.0 / 24.0);
            t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
            xWeight[1] = t1 + t0;
            xWeight[3] = t1 - t0;
            xWeight[4] = xWeight[0] + t0 + (1.0 / 2.0) * w;
            xWeight[2] = 1.0 - xWeight[0] - xWeight[1] - xWeight[3] - xWeight[4];
            break;
        case 5L:
            w2 = w * w;
            xWeight[5] = (1.0 / 120.0) * w * w2 * w2;
            w2 -= w;
            w4 = w2 * w2;
            w -= 1.0 / 2.0;
            t = w2 * (w2 - 3.0);
            xWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - xWeight[5];
            t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
            t1 = (-1.0 / 12.0) * w * (t + 4.0);
            xWeight[2] = t0 + t1;
            xWeight[3] = t0 - t1;
            t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
            t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
            xWeight[1] = t0 + t1;
            xWeight[4] = t0 - t1;
            break;
        case 6L:
            xWeight[0] = 1.0 / 2.0 - w;
            xWeight[0] *= xWeight[0] * xWeight[0];
            xWeight[0] *= xWeight[0] / 720.0;
            xWeight[1] = (361.0 / 192.0 - w * (59.0 / 8.0 + w
                * (-185.0 / 16.0 + w * (25.0 / 3.0 + w * (-5.0 / 2.0 + w)
                * (1.0 / 2.0 + w))))) / 120.0;
            xWeight[2] = (10543.0 / 960.0 + w * (-289.0 / 16.0 + w
                * (79.0 / 16.0 + w * (43.0 / 6.0 + w * (-17.0 / 4.0 + w
                * (-1.0 + w)))))) / 48.0;
            w2 = w * w;
            xWeight[3] = (5887.0 / 320.0 - w2 * (231.0 / 16.0 - w2
                * (21.0 / 4.0 - w2))) / 36.0;
            xWeight[4] = (10543.0 / 960.0 + w * (289.0 / 16.0 + w
                * (79.0 / 16.0 + w * (-43.0 / 6.0 + w * (-17.0 / 4.0 + w
                * (1.0 + w)))))) / 48.0;
            xWeight[6] = 1.0 / 2.0 + w;
            xWeight[6] *= xWeight[6] * xWeight[6];
            xWeight[6] *= xWeight[6] / 720.0;
            xWeight[5] = 1.0 - xWeight[0] - xWeight[1] - xWeight[2] - xWeight[3]
                - xWeight[4] - xWeight[6];
            break;
        case 7L:
            xWeight[0] = 1.0 - w;
            xWeight[0] *= xWeight[0];
            xWeight[0] *= xWeight[0] * xWeight[0];
            xWeight[0] *= (1.0 - w) / 5040.0;
            w2 = w * w;
            xWeight[1] = (120.0 / 7.0 + w * (-56.0 + w * (72.0 + w
                * (-40.0 + w2 * (12.0 + w * (-6.0 + w)))))) / 720.0;
            xWeight[2] = (397.0 / 7.0 - w * (245.0 / 3.0 + w * (-15.0 + w
                * (-95.0 / 3.0 + w * (15.0 + w * (5.0 + w
                * (-5.0 + w))))))) / 240.0;
            xWeight[3] = (2416.0 / 35.0 + w2 * (-48.0 + w2 * (16.0 + w2
                * (-4.0 + w)))) / 144.0;
            xWeight[4] = (1191.0 / 35.0 - w * (-49.0 + w * (-9.0 + w
                * (19.0 + w * (-3.0 + w) * (-3.0 + w2))))) / 144.0;
            xWeight[5] = (40.0 / 7.0 + w * (56.0 / 3.0 + w * (24.0 + w
                * (40.0 / 3.0 + w2 * (-4.0 + w * (-2.0 + w)))))) / 240.0;
            xWeight[7] = w2;
            xWeight[7] *= xWeight[7] * xWeight[7];
            xWeight[7] *= w / 5040.0;
            xWeight[6] = 1.0 - xWeight[0] - xWeight[1] - xWeight[2] - xWeight[3]
                - xWeight[4] - xWeight[5] - xWeight[7];
            break;
        case 8L:
            xWeight[0] = 1.0 / 2.0 - w;
            xWeight[0] *= xWeight[0];
            xWeight[0] *= xWeight[0];
            xWeight[0] *= xWeight[0] / 40320.0;
            w2 = w * w;
            xWeight[1] = (39.0 / 16.0 - w * (6.0 + w * (-9.0 / 2.0 + w2)))
                * (21.0 / 16.0 + w * (-15.0 / 4.0 + w * (9.0 / 2.0 + w
                * (-3.0 + w)))) / 5040.0;
            xWeight[2] = (82903.0 / 1792.0 + w * (-4177.0 / 32.0 + w
                * (2275.0 / 16.0 + w * (-487.0 / 8.0 + w * (-85.0 / 8.0 + w
                * (41.0 / 2.0 + w * (-5.0 + w * (-2.0 + w)))))))) / 1440.0;
            xWeight[3] = (310661.0 / 1792.0 - w * (14219.0 / 64.0 + w
                * (-199.0 / 8.0 + w * (-1327.0 / 16.0 + w * (245.0 / 8.0 + w
                * (53.0 / 4.0 + w * (-8.0 + w * (-1.0 + w)))))))) / 720.0;
            xWeight[4] = (2337507.0 / 8960.0 + w2 * (-2601.0 / 16.0 + w2
                * (387.0 / 8.0 + w2 * (-9.0 + w2)))) / 576.0;
            xWeight[5] = (310661.0 / 1792.0 - w * (-14219.0 / 64.0 + w
                * (-199.0 / 8.0 + w * (1327.0 / 16.0 + w * (245.0 / 8.0 + w
                * (-53.0 / 4.0 + w * (-8.0 + w * (1.0 + w)))))))) / 720.0;
            xWeight[7] = (39.0 / 16.0 - w * (-6.0 + w * (-9.0 / 2.0 + w2)))
                * (21.0 / 16.0 + w * (15.0 / 4.0 + w * (9.0 / 2.0 + w
                * (3.0 + w)))) / 5040.0;
            xWeight[8] = 1.0 / 2.0 + w;
            xWeight[8] *= xWeight[8];
            xWeight[8] *= xWeight[8];
            xWeight[8] *= xWeight[8] / 40320.0;
            xWeight[6] = 1.0 - xWeight[0] - xWeight[1] - xWeight[2] - xWeight[3]
                - xWeight[4] - xWeight[5] - xWeight[7] - xWeight[8];
            break;
        case 9L:
            xWeight[0] = 1.0 - w;
            xWeight[0] *= xWeight[0];
            xWeight[0] *= xWeight[0];
            xWeight[0] *= xWeight[0] * (1.0 - w) / 362880.0;
            xWeight[1] = (502.0 / 9.0 + w * (-246.0 + w * (472.0 + w
                * (-504.0 + w * (308.0 + w * (-84.0 + w * (-56.0 / 3.0 + w
                * (24.0 + w * (-8.0 + w))))))))) / 40320.0;
            xWeight[2] = (3652.0 / 9.0 - w * (2023.0 / 2.0 + w * (-952.0 + w
                * (938.0 / 3.0 + w * (112.0 + w * (-119.0 + w * (56.0 / 3.0 + w
                * (14.0 + w * (-7.0 + w))))))))) / 10080.0;
            xWeight[3] = (44117.0 / 42.0 + w * (-2427.0 / 2.0 + w * (66.0 + w
                * (434.0 + w * (-129.0 + w * (-69.0 + w * (34.0 + w * (6.0 + w
                * (-6.0 + w))))))))) / 4320.0;
            w2 = w * w;
            xWeight[4] = (78095.0 / 63.0 - w2 * (700.0 + w2 * (-190.0 + w2
                * (100.0 / 3.0 + w2 * (-5.0 + w))))) / 2880.0;
            xWeight[5] = (44117.0 / 63.0 + w * (809.0 + w * (44.0 + w
                * (-868.0 / 3.0 + w * (-86.0 + w * (46.0 + w * (68.0 / 3.0 + w
                * (-4.0 + w * (-4.0 + w))))))))) / 2880.0;
            xWeight[6] = (3652.0 / 21.0 - w * (-867.0 / 2.0 + w * (-408.0 + w
                * (-134.0 + w * (48.0 + w * (51.0 + w * (-4.0 + w) * (-1.0 + w)
                * (2.0 + w))))))) / 4320.0;
            xWeight[7] = (251.0 / 18.0 + w * (123.0 / 2.0 + w * (118.0 + w
                * (126.0 + w * (77.0 + w * (21.0 + w * (-14.0 / 3.0 + w
                * (-6.0 + w * (-2.0 + w))))))))) / 10080.0;
            xWeight[9] = w2 * w2;
            xWeight[9] *= xWeight[9] * w / 362880.0;
            xWeight[8] = 1.0 - xWeight[0] - xWeight[1] - xWeight[2] - xWeight[3]
                - xWeight[4] - xWeight[5] - xWeight[6] - xWeight[7] - xWeight[9];
            break;
        default:
            return(0);
    }
    return(1);
} /* end GetInterpolationWeights */

/*--------------------------------------------------------------------------*/
template<class T>
int vtkImageBSplineCoefficients::
InterpolatedValue
    (
        const T *Bcoeff,      /* input B-spline array of coefficients */
        T       *v,           /* resulting pixel value */
        long    Width,        /* width of the image */
        long    Height,       /* height of the image */
        long    Slices,       /* number of slices of the image */
        long    Depth,        /* number of samples per pixel */
        double  x,            /* x coordinate where to interpolate */
        double  y,            /* y coordinate where to interpolate */
        double  z,            /* y coordinate where to interpolate */
        long    SplineDegree, /* degree of the spline model */
        long    Border        /* what to do at the border */
    )

{ /* begin InterpolatedValue */

    const T *p1, *p2, *p3;
    double  xWeight[10], yWeight[10], zWeight[10];
    double  interpolated;
    double  w, u;
    double  s, t, r;
    long    xIndex[10], yIndex[10], zIndex[10];
    ptrdiff_t xIncrement, yIncrement, zIncrement;
    long    Width2 = 2L * Width - 2L;
    long    Height2 = 2L * Height - 2L;
    long    Slices2 = 2L * Slices - 2L;
    long    CentralIndex = SplineDegree / 2L;
    long    i, j, k, l, c;
    long    imax, jmax, kmax;

    if (SplineDegree < 2 || SplineDegree > 9)
      {
      return 0;
      }

    /* check for 1D and 2D images */
    imax = SplineDegree;
    jmax = SplineDegree;
    kmax = SplineDegree;
    if (Width == 1)
      {
      imax = 0;
      }
    if (Height == 1)
      {
      jmax = 0;
      }
    if (Slices == 1)
      {
      kmax = 0;
      }

    /* compute the interpolation indexes and fractional offsets*/
    if (SplineDegree & 1L) {
        i = (long)(floor(x));
        j = (long)(floor(y));
        k = (long)(floor(z));
    }
    else {
        i = (long)(floor(x + 0.5));
        j = (long)(floor(y + 0.5));
        k = (long)(floor(z + 0.5));
    }

    s = x - i;
    t = y - j;
    r = z - k;

    i -= CentralIndex;
    j -= CentralIndex;
    k -= CentralIndex;

    for (l = 0L; l <= SplineDegree; l++) {
        xIndex[l] = i++;
        yIndex[l] = j++;
        zIndex[l] = k++;
    }

    /* get the interpolation weights */
    xWeight[0] = 1.0;
    yWeight[0] = 1.0;
    zWeight[0] = 1.0;

    if (Width > 1) {
        GetInterpolationWeights(xWeight, s, SplineDegree);
    }
    if (Height > 1) {
        GetInterpolationWeights(yWeight, t, SplineDegree);
    }
    if (Slices > 1) {
        GetInterpolationWeights(zWeight, r, SplineDegree);
    }

    switch (Border)
        {
        case VTK_BSPLINE_CLAMP:
            /* apply the constant boundary conditions */
            for (l = 0L; l <= SplineDegree; l++) {
                if (xIndex[l] < 0) {
                    xIndex[l] = 0;
                }
                else if (xIndex[l] >= Width) {
                    xIndex[l] = Width - 1;
                }
            }
            for (l = 0L; l <= SplineDegree; l++) {
                if (yIndex[l] < 0) {
                    yIndex[l] = 0;
                }
                else if (yIndex[l] >= Height) {
                    yIndex[l] = Height - 1;
                }
            }
            for (l = 0L; l <= SplineDegree; l++) {
                if (zIndex[l] < 0) {
                    zIndex[l] = 0;
                }
                else if (zIndex[l] >= Slices) {
                    zIndex[l] = Slices - 1;
                }
            }
            break;

        case VTK_BSPLINE_MIRROR:
            /* apply the mirror boundary conditions */
            for (l = 0L; l <= SplineDegree; l++) {
                xIndex[l] = (Width == 1L) ? (0L) : ((xIndex[l] < 0L) ?
                    (-xIndex[l] - Width2 * ((-xIndex[l]) / Width2))
                    : (xIndex[l] - Width2 * (xIndex[l] / Width2)));
                if (Width <= xIndex[l]) {
                    xIndex[l] = Width2 - xIndex[l];
                }
                yIndex[l] = (Height == 1L) ? (0L) : ((yIndex[l] < 0L) ?
                    (-yIndex[l] - Height2 * ((-yIndex[l]) / Height2))
                    : (yIndex[l] - Height2 * (yIndex[l] / Height2)));
                if (Height <= yIndex[l]) {
                    yIndex[l] = Height2 - yIndex[l];
                }
                zIndex[l] = (Slices == 1L) ? (0L) : ((zIndex[l] < 0L) ?
                    (-zIndex[l] - Slices2 * ((-zIndex[l]) / Slices2))
                    : (zIndex[l] - Slices2 * (zIndex[l] / Slices2)));
                if (Slices <= zIndex[l]) {
                    zIndex[l] = Slices2 - zIndex[l];
                }
            }
            break;

        case VTK_BSPLINE_REPEAT:
            /* apply the repeat boundary conditions */
            for (l = 0L; l <= SplineDegree; l++) {
                if ((xIndex[l] %= Width) < 0) {
                    xIndex[l] += Width;
                }
                if ((yIndex[l] %= Height) < 0) {
                    yIndex[l] += Height;
                }
                if ((zIndex[l] %= Slices) < 0) {
                    zIndex[l] += Slices;
                }
            }
            break;
    }

    /* precompute the increments in each direction */
    xIncrement = Depth;
    yIncrement = xIncrement * Width;
    zIncrement = yIncrement * Height;

    /* perform interpolation */
    for (c = 0L; c < Depth; c++) {
        p1 = Bcoeff + c;
        interpolated = 0.0;
        for (k = 0L; k <= kmax; k++) {
            p2 = p1 + (zIndex[k] * zIncrement);
            u = 0.0;
            for (j = 0L; j <= jmax; j++) {
                p3 = p2 + (yIndex[j] * yIncrement);
                w = 0.0;
                for (i = 0L; i <= imax; i++) {
                    w += xWeight[i] * p3[xIndex[i] * xIncrement];
                }
                u += yWeight[j] * w;
            }
            interpolated += zWeight[k] * u;
        }
        v[c] = interpolated;
    }

    return(1);
} /* end InterpolateValue */
