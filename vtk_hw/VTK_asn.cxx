#include "vtkSmartPointer.h"
#include "vtkImageSliceMapper.h"
#include "vtkImageActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkCamera.h"
#include "vtkInteractorStyleImage.h"
#include "vtkCommand.h"
#include "vtkImageSlice.h"
#include "vtkImageProperty.h"
#include "vtkDiscreteMarchingCubes.h"
#include "vtkPolyDataMapper.h"
#include "vtkPoints.h"
#include "vtkNew.h"
#include "vtkProgrammableFilter.h"
#include "vtkCallbackCommand.h"
#include "vtkProperty.h"
#include "vtkSphereSource.h"
#include "vtkPolyPointSource.h"
#include "vtkTextActor.h"
#include "vtkTextProperty.h"
#include "vtkPolyDataNormals.h"
#include "vtkArrowSource.h"
#include "vtkGlyph3D.h"
#include "vtkDistancePolyDataFilter.h"
#include "vtkCommand.h"
#include "vtkRendererCollection.h"
#include "vtkCellPicker.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageToVTKImageFilter.h"
#include "itkFlipImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkDiffusionTensor3D.h"
#include "itkSymmetricEigenAnalysisImageFilter.h"
#include "itkTensorFractionalAnisotropyImageFilter.h"
#include "itkVector.h"
#include "itkImageIterator.h"

const unsigned int nDims = 3 ;   
typedef itk::DiffusionTensor3D < double > DiffTensor ;
typedef itk::Image < DiffTensor , nDims > ImageType ;
typedef itk::Vector <double, nDims> VectorType;
typedef itk::Image <VectorType, nDims> PAImageType ;
typedef itk::Image < double, nDims> FAImageType ;
typedef itk::Image < int, nDims> TAImg ;
typedef itk::Image < int, nDims> SegImg ;

vtkNew<vtkPoints> points;
std::list<ImageType::IndexType> img_ls;


struct track_ls
{
  std::list<ImageType::IndexType> trackSeedls;
  std::list<ImageType::IndexType> segmenttrackPixelList;
} ;

vtkSmartPointer<vtkPolyData> DecimateMesh(vtkSmartPointer<vtkPolyData> inputMesh, double ratio)
{
  vtkSmartPointer<vtkDecimatePro> decimateFilter = vtkSmartPointer<vtkDecimatePro>::New();
  decimateFilter->SetInputData(inputMesh);
  decimateFilter->SetTargetReduction(ratio);
  decimateFilter->Update();

  return decimateFilter->GetOutput();
}

ImageType::IndexType getCurLoc(VectorType curVec, double delta, ImageType::IndexType currLoc, bool sign){
  ImageType::IndexType newLoc = currLoc;
  int sign = 1;
  if (!sign){
    sign = -1
  }
  newLoc[0] = round(sign*curVec[0]*delta + currLoc[0]);
  newLoc[1] = round(sign*curVec[1]*delta + currLoc[1]);
  newLoc[2] = round(sign*curVec[2]*delta + currLoc[2]);
  
  return newLoc;
}


int trackRe(PAImageType::Pointer getPAImg, FAImageType::Pointer myFAImage, TAImg::Pointer myTrackImage, ImageType::IndexType thisloc, double delta, int iter, std::list<ImageType::IndexType> & trackPixelList)
{
  if (!myTrackImage->GetLargestPossibleRegion().IsInside(thisloc)) return 0;
  if (myTrackImage->GetPixel(thisloc) == 1.0) return 0 ;
  if (myFAImage->GetPixel(thisloc) < 0.25) return 0;
  if (iter > 1000) return 0;
  myTrackImage->SetPixel(thisloc , 1.0);
  trackPixelList.push_back(thisloc);
  iter++;
  ImageType::IndexType newlocplus = thisloc;
  ImageType::IndexType newlocneg = thisloc;
  VectorType curVec;
  curVec = getPAImg->GetPixel(thisloc);
  ImageType::IndexType newLoc = getNewIdx(curVec, delta, thisLoc, true);
  ImageType::IndexType newLocNeg = getNewIdx(curVec, delta, thisLoc, false);
  trackRe(getPAImg, myFAImage, myTrackImage, newLoc, delta, iter, trackPixelList);
  trackRe(getPAImg, myFAImage, myTrackImage, newLocNeg , delta, iter, trackPixelList);

  return 0;
}

FAImageType::Pointer getFAImg(char * fileCheck) 
{
  typedef itk::ImageFileReader < ImageType > ImageReaderType1 ;
  ImageReaderType1 ::Pointer myReader1 = ImageReaderType1 ::New() ;  
  myReader1->SetFileName (fileCheck) ;
  myReader1->Update();
  ImageType::Pointer Imge = myReader1->GetOutput();
  
  typedef itk::TensorFractionalAnisotropyImageFilter <ImageType, FAImageType> FAFilterType ;
  FAFilterType::Pointer myFAImageFilter = FAFilterType::New() ;
  myFAImageFilter ->SetInput(Imge) ;
  myFAImageFilter ->Update() ;

  return myFAImageFilter->GetOutput();
}

void mouseClick()
{
  int* pos = this->GetInteractor()->GetEventPosition();

  vtkSmartPointer<vtkCellPicker> picker =
    vtkSmartPointer<vtkCellPicker>::New();
  picker->SetTolerance(/*0.0005*/ 0.001);
  picker->Pick(pos[0], pos[1], 0, this->GetDefaultRenderer());

  double* worldPosition = picker->GetPickPosition();
  std::cout << "check " << picker->GetCellId() << std::endl;

  if(picker->GetCellId() != -1)
  {
    std::cout << "position " << worldPosition[0] << " " << worldPosition[1]
    << " " << worldPosition[2] << endl;
    vtkSmartPointer<vtkIdTypeArray> ids = vtkSmartPointer<vtkIdTypeArray>::New();
    ids->SetNumberOfComponents(1);
    ids->InsertNextValue(picker->GetCellId());
    vtkSmartPointer<vtkSelectionNode> selectionNode = vtkSmartPointer<vtkSelectionNode>::New();
    selectionNode->SetFieldType(vtkSelectionNode::CELL);
    selectionNode->SetContentType(vtkSelectionNode::INDICES);
    selectionNode->SetSelectionList(ids);
    vtkSmartPointer<vtkSelection> selection = vtkSmartPointer<vtkSelection>::New();
    selection->AddNode(selectionNode);
  }
}

namespace
{
  void folTimer(vtkObject *caller, long unsigned int eventId,
                             void *clientData, void *callData);

  unsigned int counter = 0;
  unsigned int maxCount = 700000;

  void checkDot(void *arguments)
  {
    vtkProgrammableFilter *filter = static_cast<vtkProgrammableFilter *>(arguments);

    vtkPoints *pointIn = filter->GetPolyDataInput()->GetPoints();
    vtkIdType numPts = pointIn->GetNumberOfPoints();
    vtkNew<vtkPoints> newPts;
    newPts->SetNumberOfPoints(numPts);

    double p[3];
    for (int i = 0; i < numPts; i++)
    {
      if (i < counter + 1)
      {
        points->GetPoint(i, p);
        newPts->SetPoint(i, p);
      }
      else
      {
        pointIn->GetPoint(i, p);
        newPts->SetPoint(i, p);
      }
    }

    filter->GetPolyDataOutput()->CopyStructure(filter->GetPolyDataInput());
    filter->GetPolyDataOutput()->SetPoints(newPts);
    std::cout << "Iteration: " << counter << ", Number of Points: " << numPts << std::endl;
  }

} 

namespace
{
  void folTimer(vtkObject *caller, long unsigned int vtkNotUsed(eventId), void *clientData, void *vtkNotUsed(callData))
  {

    auto filter = static_cast<vtkProgrammableFilter *>(clientData);

    vtkRenderWindowInteractor *renderO =
        static_cast<vtkRenderWindowInteractor *>(caller);

    filter->Modified();

    renderO->Render();

    if (counter > maxCount)
      renderO->DestroyTimer();

    counter++;
  }
}

TAImg::Pointer setImg(ImageType::Pointer Imge, ImageType::RegionType myReigon){
  TAImg::Pointer myTrackImage = TAImg::New() ;
  myTrackImage->SetOrigin(Imge->GetOrigin() ) ;
  myTrackImage->SetDirection(Imge->GetDirection() );
  myTrackImage->SetSpacing(Imge->GetSpacing() );
  myTrackImage->SetRegions(myRegion);
  return myTrackImage;
}

track_ls setTracker(char * fileCheck,  char * SegmentCCFileName, ImageType::IndexType seeds)
{
  typedef itk::ImageFileReader < ImageType > ImageReaderType1 ;
  ImageReaderType1 ::Pointer myReader1 = ImageReaderType1 ::New() ;  
  myReader1->SetFileName (fileCheck) ;
  myReader1->Update();   // go read
  ImageType::Pointer Imge = myReader1->GetOutput();

  typedef itk::ImageFileReader < SegImg > ImageReaderType2 ;
  ImageReaderType2 ::Pointer myReader2 = ImageReaderType2::New() ;  
  myReader2 ->SetFileName ( SegmentCCFileName ) ;
  myReader2 ->Update();   // go read
  SegImg::Pointer SegmentedCCImage = myReader2->GetOutput();
   
  ImageType::RegionType regionCus;
  ImageType::SizeType size = Imge->GetLargestPossibleRegion().GetSize() ;
  ImageType::IndexType idx ;
  idx[0] = 0; 
  idx[1] = 0; 
  idx[2] = 0;
  
  regionCus.SetSize( size ) ;
  regionCus.SetIndex( idx ) ;
  typedef itk::ImageRegionIterator < ImageType > InputIteratorType ;
  typedef itk::ImageRegionIterator < PAImageType > OutputIteratorType ;

  PAImageType::Pointer getPAImg = setImg(Imge, regionCus);
  getPAImg->Allocate() ;

  OutputIteratorType outputIterator (getPAImg, regionCus);
  InputIteratorType inputIterator (Imge, regionCus);
  outputIterator.GoToBegin ();
  inputIterator.GoToBegin () ;
  DiffTensor thisTensor ;
  DiffTensor::EigenValuesArrayType myEVAT;
  DiffTensor::EigenVectorsMatrixType myEVmatr ;
  VectorType curVec ;

  while (!outputIterator.IsAtEnd() )
  {
   thisTensor =  inputIterator.Value() ;
   thisTensor.ComputeEigenAnalysis(myEVAT, myEVmatr) ;
   curVec[0] = myEVmatr[2][0]*1 ; curVec[1] = myEVmatr[2][1]*1 ; curVec[2] = myEVmatr[2][2]*1 ; // Principal axis vector

   if (myEVmatr[2][2] == 1) { 
     curVec[0] = 0; curVec[1] = 0; curVec[2] = 0; //Change zero tensor to 0 Principal Vector Direction
   }
   outputIterator.Set(curVec) ;
   ++outputIterator ;
   ++inputIterator ;
  }

  typedef itk::TensorFractionalAnisotropyImageFilter <ImageType, FAImageType> FAFilterType ;
  FAFilterType::Pointer myFAImageFilter = FAFilterType::New() ;
  myFAImageFilter ->SetInput(Imge) ;
  myFAImageFilter ->Update() ;

  TAImg::Pointer myTrackImage = setImg(Imge, regionCus);
  myTrackImage->Allocate() ;

  int iter = 0;
  double delta = 1.5;
  ImageType::IndexType thisloc = seeds;
  ImageType::IndexType newloc = seeds;
  
  std::list<ImageType::IndexType> trackSeedls;
  trackRe(getPAImg, myFAImageFilter->GetOutput(), myTrackImage, thisloc, delta, iter, trackSeedls);
  
  TAImg::Pointer myCCTrackImage = setImg(Imge, regionCus);
  myCCTrackImage ->Allocate() ;

  typedef itk::ImageRegionIterator < SegImg > SegmentIteratorType ;
  SegmentIteratorType segmentIterator (SegmentedCCImage, regionCus);
  segmentIterator .GoToBegin ();
  
  std::list<ImageType::IndexType> segmenttrackPixelList;

  while (!segmentIterator.IsAtEnd() )
  {
   if (segmentIterator.Value() == 1.0)
   {
      int iter = 0;
      trackRe(getPAImg, myFAImageFilter->GetOutput(), myCCTrackImage , segmentIterator.GetIndex(), delta, iter, segmenttrackPixelList);
   }
   ++segmentIterator ;
  }

  track_ls result = {trackSeedls, segmenttrackPixelList};
  return result;
}


 TAImg::Pointer buildITK(FAImageType::Pointer myITKImage, ){
  TAImg::Pointer myTrackImage = setImg(myITKImage, myRegion);
  myTrackImage->Allocate() ;
  myTrackImage->SetPixel(seeds , 1.0);
  return myTrackImage;
 }

 vtkSmartPointer < vtkImageSliceMapper > buildMapper(filter){
  filterImg::Pointer vtkSmartPointer < vtkImageSliceMapper > imgMap = vtkSmartPointer < vtkImageSliceMapper > ::New() ;
  imgMap->SetInputData ( filter->GetOutput() ) ;
  imgMap->SetOrientationToX () ;
  imgMap->SetSliceNumber ( 55 ) ;
  imgMap->SliceAtFocalPointOn () ;
  imgMap->SliceFacesCameraOn () ;  
  return imgMap;
 }



int main ( int argc, char * argv[] )
{
  // Verify command line arguments
  if( argc < 4 )
    {
      std::cerr << "Usage: " << std::endl ;
      std::cerr << argv[0] << " FAImg " << std::endl ; 
      std::cerr <<" seeds, segmentation track" << std::endl;
      return -1 ;
    }
  int seedOrSeg = strtol(argv[3], NULL, 10);
  FAImageType::Pointer myITKImage;
  typedef itk::ImageFileReader < ImageType > ImageReaderType1 ;
  ImageReaderType1 ::Pointer myReader1 = ImageReaderType1 ::New() ;  
  myReader1->SetFileName (argv[1]) ;
  myReader1->Update();
  ImageType::Pointer Imge = myReader1->GetOutput();
  
  typedef itk::TensorFractionalAnisotropyImageFilter <ImageType, FAImageType> FAFilterType ;
  FAFilterType::Pointer myFAImageFilter = FAFilterType::New() ;
  myFAImageFilter ->SetInput(Imge) ;
  myFAImageFilter ->Update() ;

  myITKImage = myFAImageFilter->GetOutput();

  ImageType::IndexType seeds ;
  seeds[0] = 75; 
  seeds[1] = 83; 
  seeds[2] = 52;
  track_ls myResults = setTracker(argv[1], argv[2], seeds); 

  ImageType::RegionType regionCus;
  ImageType::SizeType size = myITKImage->GetLargestPossibleRegion().GetSize() ;
  ImageType::IndexType idx ;
  idx[0] = 0; 
  idx[1] = 0; 
  idx[2] = 0;
  regionCus.SetSize( size ) ;
  regionCus.SetIndex( idx ) ;
  
  TAImg::Pointer myTrackImage = setImg(myITKImage, regionCus);
  myTrackImage->Allocate() ;
  myTrackImage->SetPixel(seeds , 1.0);

  typedef itk::FlipImageFilter < FAImageType > FIFType ;
  FIFType::Pointer myFlip = FIFType::New();
  myFlip->SetInput(myITKImage);
  FIFType::FlipAxesArrayType reverseAxe;
  reverseAxe[0] = true; 
  reverseAxe[1] = true; 
  reverseAxe[2] = false;
  myFlip->SetFlipAxes(reverseAxe);

  typedef itk::ImageToVTKImageFilter < FAImageType > ITKToVTKFilterType ;
  ITKToVTKFilterType::Pointer itkToVTKfilter = ITKToVTKFilterType::New() ;
  itkToVTKfilter->SetInput ( myFlip->GetOutput() ) ;
  itkToVTKfilter->Update() ;
  
  typedef itk::ImageToVTKImageFilter < SegImg > filterImg ;
  filterImg::Pointer itkMaskToVTKfilter = filterImg::New() ;
  itkMaskToVTKfilter->SetInput ( myTrackImage) ;
  itkMaskToVTKfilter->Update() ;

  vtkSmartPointer < vtkImageSliceMapper > imgMap = buildMapper(itkToVTKfilter);
  vtkSmartPointer <vtkImageSlice> slice = vtkSmartPointer <vtkImageSlice> ::New();
  slice->SetMapper(imgMap);
  slice->GetProperty()->SetColorWindow(0.8877);
  slice->GetProperty()->SetColorLevel(0.4446);

  if (seedOrSeg == 1){
      img_ls = myResults.trackSeedls;
  }
  if (seedOrSeg == 2){
      img_ls = myResults.segmenttrackPixelList;
  }
  vtkNew<vtkPolyPointSource> polypointSource;

  ImageType::IndexType thisIndex;
  for (std::list<ImageType::IndexType>::iterator it = img_ls.begin(); it!=img_ls.end(); it++){
    thisIndex = *it;
    points->InsertNextPoint(thisIndex[0], thisIndex[1], thisIndex[2]);
  }
  
  vtkNew<vtkPoints> initialPoints;
  std::list<ImageType::IndexType>::iterator it = img_ls.begin();
  thisIndex = *it;
  for (std::list<ImageType::IndexType>::iterator it = img_ls.begin(); it!=img_ls.end(); it++){
    initialPoints->InsertNextPoint(thisIndex[0], thisIndex[1], thisIndex[2]);
  }


  polypointSource->SetPoints(initialPoints);
  polypointSource->Update();
  vtkNew<vtkProgrammableFilter> filter;
  filter->SetInputConnection(polypointSource->GetOutputPort());
  filter->SetExecuteMethod(checkDot, filter);
  
  vtkSmartPointer < vtkPolyDataMapper > polyMapper = vtkSmartPointer < vtkPolyDataMapper > ::New() ;
  polyMapper->SetInputConnection (filter->GetOutputPort()) ;
  vtkSmartPointer < vtkActor > actor = vtkSmartPointer < vtkActor >::New() ;
  actor->SetMapper(polyMapper);

  vtkSmartPointer < vtkTextActor > text1 = vtkSmartPointer < vtkTextActor >::New() ;
  text1->SetInput ( "FA Image" ) ;
  text1->SetDisplayPosition ( 50, 450 ) ;
  text1->GetTextProperty()->SetColor ( 1, 1, 1 ) ;
  text1->GetTextProperty()->SetFontSize ( 20 ) ;

  vtkSmartPointer < vtkTextActor > text2 = vtkSmartPointer < vtkTextActor >::New() ;
  if (seedOrSeg == 1){
    text2->SetInput ( "seeds check" ) ;
  }
  if (seedOrSeg == 2){
      text2->SetInput ( "segmentation check" ) ;
  }
  text2->SetDisplayPosition ( 500, 500 ) ;
  text2->GetTextProperty()->SetColor ( 1, 1, 1 ) ;
  text2->GetTextProperty()->SetFontSize ( 30 ) ;

  vtkSmartPointer < vtkRenderer > renderer1 = vtkSmartPointer < vtkRenderer >::New() ;
  renderer1->AddViewProp ( slice ) ;
  renderer1->SetViewport(0,0,0.5,1);
  renderer1->AddActor(text1);

  vtkSmartPointer < vtkCamera > camera = renderer1->GetActiveCamera() ;

  double position[3],  imageCenter[3] ;
  itkToVTKfilter->GetOutput()->GetCenter(imageCenter);
  position[0] = imageCenter[0] ;
  position[1] = imageCenter[1] ;
  position[2] = -160 ;
  double spacing[3] ;
  int imageDims[3] ;
  itkToVTKfilter->GetOutput()->GetSpacing ( spacing ) ;
  itkToVTKfilter->GetOutput()->GetDimensions ( imageDims ) ;
  double imagePhysicalSize[3] ;
  for ( unsigned int d = 0 ; d < 3 ; d++ ) imagePhysicalSize[d] = spacing[d] * imageDims[d] ;
  camera->ParallelProjectionOn(); 
  camera->SetFocalPoint(imageCenter); 
  camera->SetPosition(position);
  camera->SetParallelScale(imageDims[2] / 0.8);

  vtkSmartPointer < vtkRenderer > renderer2 = vtkSmartPointer < vtkRenderer >::New() ;
  renderer2->SetViewport(0.5, 0, 1, 1);
  renderer2->AddActor(actor);
  renderer2->AddActor(text2);

  double thisIndexFP[3], thisIndexPosition[3];
  thisIndexFP[0] = thisIndex[0]; 
  thisIndexFP[1] = thisIndex[1]; 
  thisIndexFP[2] = thisIndex[2];
  thisIndexPosition[0] = thisIndexFP[0];
  thisIndexPosition[1] = thisIndexFP[1];
  thisIndexPosition[2] = -160;
  renderer2->ResetCamera();
  camera = renderer2->GetActiveCamera();
  camera->ParallelProjectionOn(); camera->SetFocalPoint(thisIndexFP); camera->SetPosition(thisIndexPosition);
  camera->SetParallelScale ( imageDims[2] / 0.8 ) ;

  vtkSmartPointer < vtkRenderWindow > window = vtkSmartPointer < vtkRenderWindow >::New() ;
  vtkSmartPointer < vtkInteractorStyleImage > style = vtkSmartPointer < vtkInteractorStyleImage >::New() ;
  vtkSmartPointer < vtkRenderWindowInteractor > interactor = vtkSmartPointer < vtkRenderWindowInteractor >::New() ;
  window->AddRenderer ( renderer1 ); window->AddRenderer ( renderer2 ); window->SetSize ( 1000, 500 ) ;
  window->Render();
  interactor->SetRenderWindow ( window ) ;
  style->SetInteractionModeToImageSlicing();
  interactor->Initialize() ;
  interactor->CreateRepeatingTimer ( 100 ) ;

  vtkNew<vtkCallbackCommand> timerCallback;
  timerCallback->SetCallback(folTimer);
  timerCallback->SetClientData(filter);

  interactor->AddObserver(vtkCommand::TimerEvent, timerCallback);

  interactor->Start() ;
  interactor->ReInitialize();

  return 0 ;
}


