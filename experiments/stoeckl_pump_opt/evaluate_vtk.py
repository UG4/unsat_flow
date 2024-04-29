"""Script to extract the lens thickness from Stoeckl Lens Testcase using Paraview's Python interface"""

import sys
sys.path.append("/opt/paraview/lib/python3.12/site-packages/")

import paraview.simple as pv
import glob

def collectFiles(pathToFiles: str) -> list:
    """Collect all vtk files to be processed"""
    return sorted(glob.glob(pathToFiles))[-2:-1]

def extractToeLength(filename: str):
    """Extract Lens Thickness from the vtk files. 

    Steps: 
        1. Read the vtk file into Paraview
        2. Apply Contour filter on concentration with value 0.5
        3. Export the contour points
    """
    toeLength = []
    upperBound = 0.3

    pv.Connect()
    vtkFiles = collectFiles("./*.pvtu")

    if len(vtkFiles) == 0:
        vtkFiles = collectFiles("./*.vtu")
    
    reader = pv.OpenDataFile(vtkFiles)
    if reader:
        print("Success")
    else:
        print("Failed")

    contour = pv.Contour(Input=reader)
    contour.ContourBy = "c"
    contour.Isosurfaces = [0.5]

    pv.UpdatePipeline()

    animationScene = pv.GetAnimationScene()
    animationScene.GoToLast()

    pv.UpdatePipeline(proxy=contour)
    data = pv.servermanager.Fetch(contour)
    
    if data.GetPoints() is None:
        #print("No data found for timestep ", t)
        pass
    
    bounds = data.GetPoints().GetBounds()
    # bounds are [xmin, xmax, ymin, ymax, zmin, zmax]

    with open(f"{filename}", "w") as f:
        f.write("time,value\n")
        f.write(f"{35000},{bounds[1]}\n")
        f.write(f"{35000},{bounds[1]}\n")
        f.write("FINISHED,FINISHED")

if __name__ == "__main__":
    filename = sys.argv[1]

    extractToeLength(filename)
    print("Exported Toe Length to ", filename)
