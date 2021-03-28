# Mapping-FDCs-in-GCs
Macro utilized in J.C. Anania et al 2020 Frontiers in Immunology A novel image analysis approach reveals a role for complement receptors 1 and 2 in follicular dendritic cell (FDC) organization in germinal centers (GCs)

Image J macro defines splenic white pulp, B cell follicle, GCs, LZ and DZ regions of interest (ROI) and maps distribution of FDC-M1 poditive FDCs with ROIs. FDC distance from ROI centre/edge (EDM) and organization (cluster analysis) were measured.


// FDC GC ZONE MACRO
 
run("Close All"); 
print("\\Close"); 
roiManager("Reset"); 

//***SETUP***

//Open and split channels 

//File location= C:/Users/jessa/Desktop/GC Zone Macro/
//Save location= C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/

open("C:/Users/jessa/Desktop/GC Zone Macro/SRBC11_WT_d6_1-MOMA-PE, FDCM1-BV421, B220-FITC, Ki67-e660_1.lsm");

getPixelSize(unit, pixelWidth, pixelHeight); //removing um and only working in pixels from here on
OriginalScale= pixelWidth;
run("Set Scale...", "distance=1 known=1 unit=pixels");
run("Conversions...", " ");
 
ImageIDoriginal=getImageID(); 
print(ImageIDoriginal); 
 
run("Duplicate...", "duplicate"); 
ImageIDduplicate=getImageID();
print(ImageIDduplicate);

//Split channels and rename 

run("Split Channels");
 
MOMA=ImageIDduplicate-2;//assign channels 
print(MOMA); 
selectImage(MOMA); 
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/MOMA PE.tif"); 
 
Ki67=ImageIDduplicate-1; 
print(Ki67);
selectImage(Ki67);
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/Ki67 660.tif");

B220=ImageIDduplicate-3; 
print(B220);
selectImage(B220);
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/B220 FITC.tif");
 
FDCM1=ImageIDduplicate-4;
print(FDCM1);
selectImage(FDCM1);
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/FDCM1 421.tif");
run("Duplicate...", "FDCduplicate");
rename("FDCduplicate");

//Threshold smooth and create mask
 
selectImage(MOMA); //smooth and mask 
run("Options...", "iterations=4 count=1 black do=Nothing"); 
run("Smooth"); 
setAutoThreshold("Otsu dark"); 
run("Convert to Mask"); 
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/MOMA mask.tif");  //Save for MZ mask

run("Duplicate...", "FDCduplicate"); // to generate one for the total GC
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/WP mask.tif"); 

selectImage(Ki67); 
run("Smooth"); 
setAutoThreshold("Otsu dark"); 
run("Convert to Mask"); 
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/Ki67 mask.tif");

selectImage(FDCM1);
run("Smooth");  
setAutoThreshold("RenyiEntropy dark");
run("Convert to Mask");
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/FDC mask.tif");

selectImage(B220);
run("Options...", "iterations=4 count=1 black do=Nothing"); 
run("Smooth"); 
setAutoThreshold("Triangle dark"); 
run("Convert to Mask"); 
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/B220 mask.tif");
run("Duplicate...", "duplicate"); 
B220IDduplicate=getImageID();

//***WHITE PULP***

//Identify white pulp region (convex hull)
 
selectWindow("WP mask.tif");
run("8-bit"); 
run("Options...", "iterations=3 count=1 black edm=8-bit do=Close"); 
run("Fill Holes"); 
run("Options...", "iterations=3 count=1 black edm=8-bit do=Dilate"); 
run("Analyze Particles...", "size=100-Infinity display summarize add"); 
roiManager("Select", 0); 
roiManager("Rename", "White pulp");

selectWindow("WP mask.tif");
setMinAndMax(0, 255);
run("Select None");
run("Analyze Particles...", "size=100-Infinity show=[Count Masks] display summarize add");
run("glasbey");
run("Enhance Contrast", "saturated=0.35");

run("Point Tool...");
setTool("multipoint");
waitForUser("Select MZ") // Select delete in secondary window
getSelectionCoordinates(xpoints, ypoints)

for (o=0; o<xpoints.length; o=o+1){
	value=getPixel(xpoints[o], ypoints[o]);
	changeValues(value, value, 255);
	print(o, value);
}

changeValues(0, 254, 0);
changeValues(256,65000,0);
run("Select None");
run("8-bit");
run("Options...", "iterations=1 count=1 black edm=8-bit do=Nothing");
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/MZ mask.tif");
run("Duplicate...", "mzduplicate");
rename("mzduplicate");
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/WP mask-1.tif");

selectWindow("MZ mask.tif");
run("8-bit");
run("Select None");
run("Grays");
setAutoThreshold("Huang dark");
run("Convert to Mask");
//run("Invert");
run("Create Selection");
roiManager("Add");
MZ= roiManager("count")-1;
roiManager("Select", MZ);
roiManager("Rename", "MZ"); 
run("Make Binary");
run("Make Inverse");
run("Convex Hull");
run("Create Mask");
run("Create Selection");
roiManager("Add");
line= roiManager("count")-1;
roiManager("Select", line);
roiManager("Rename", "WP outline"); 
run("Clear Outside");
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/WP mask.tif");

//Identify centre

Roi.getContainedPoints(xpoints, ypoints);
Array.getStatistics(xpoints, min, max, meanx, stdDev);
Array.getStatistics(ypoints, min, max, meany, stdDev);
print(meanx, meany);

//Create EDM for heat map (centre)

getDimensions(width, height, channel, slices, frames); 
newImage("conditional dilation", "8-bit black", width, height, 1);
CondilID=getImageID();
setPixel(meanx, meany, 255); 
run("Options...", "iterations=3 count=1 edm=16-bit do=Nothing");
run("Distance Map");
run("35 step");
outline= roiManager("count")-1;
roiManager("Select", outline);
run("Make Inverse");
changeValues(0, 65535, 0);
run("Make Inverse");
getRawStatistics(nPixels, mean, min, max, std, histogram); 
setMinAndMax(min, max);
print(min,max);
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/WP heatmap centre.tif") 
print("\\Clear");
getRawStatistics(nPixels, mean, min, max, std, histogramArray); // random distribution control
for (n = 0; n < histogramArray.length; n++) {
	print(n,",",histogramArray [n]);
}
print(nPixels);

waitForUser("Save values", "1. Save white pulp centre random distribution control values for graphpad analysis");

//MZ centre map

selectWindow("WP mask-1.tif");
imageCalculator("AND create 32-bit", "WP mask-1.tif","WP heatmap centre.tif"); 
selectWindow("Result of WP mask-1.tif"); 
run("35 step"); 
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/MZ heatmap.tif"); 

//FDC centre map

selectWindow("FDC mask.tif");
imageCalculator("AND create 32-bit", "FDC mask.tif","WP heatmap centre.tif"); 
selectWindow("Result of FDC mask.tif"); 
run("35 step"); 
roiManager("Select", outline);
run("Clear Outside");
run("Make Inverse");  
changeValues(0, 255, 0);
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/FDC heatmap.tif"); 

//Ki67 centre map

selectWindow("Ki67 mask.tif");
imageCalculator("AND create 32-bit", "Ki67 mask.tif","WP heatmap centre.tif"); 
selectWindow("Result of Ki67 mask.tif"); 
run("35 step"); 
roiManager("Select", outline);
run("Clear Outside");
run("Make Inverse");  
changeValues(0, 255, 0);
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/Ki67 heatmap.tif"); 

//Create EDM for heat map (distance)

newImage("DisEdge", "8-bit black", 1024, 1024, 1);
roiManager("Select", outline);
run("Clear Outside");
run("Make Inverse"); 
changeValues(0, 254, 255);
run("Distance Map");
run("35 step");
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/WP DisEdge heat map.tif");

print("\\Clear");
getRawStatistics(nPixels, mean, min, max, std, histogramArray); // random distribution control
for (n = 0; n < histogramArray.length; n++) {
	print(n,",",histogramArray [n]);
}
print(nPixels);

waitForUser("Save values", "2. Save white pulp edge random distribution control values for graphpad analysis");

//MZ edge map

selectWindow("WP mask-1.tif");
imageCalculator("AND create 32-bit", "WP mask-1.tif","WP DisEdge heat map.tif"); 
selectWindow("Result of WP mask-1.tif"); 
run("35 step"); 
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/MZ DisEdge heat map.tif"); 

//FDC edge map

selectWindow("FDC mask.tif");
roiManager("Select", outline);
run("Clear Outside");
run("Make Inverse");
changeValues(0, 255, 0);
imageCalculator("AND create 32-bit", "FDC mask.tif","WP DisEdge heat map.tif"); 
selectWindow("Result of FDC mask.tif"); 
run("35 step"); 
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/FDC WP DisEdge heat map.tif"); 

//Ki67 edge map

selectWindow("Ki67 mask.tif");
imageCalculator("AND create 32-bit", "Ki67 mask.tif","WP DisEdge heat map.tif"); 
selectWindow("Result of Ki67 mask.tif"); 
run("35 step"); 
roiManager("Select", outline);
run("Clear Outside");
run("Make Inverse");  
changeValues(0, 255, 0);
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/Ki67 WP DisEdge heatmap.tif"); 

//measure histo pixel distances

open("C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/FDCM1 421.tif");

selectWindow("FDC WP DisEdge heat map.tif");
getRawStatistics(nPixels, mean, min, max, std, dhistogram); 
selectWindow("FDC mask.tif");
selectWindow("Results")
print("nResults", nResults);
print("\\Clear");

setBatchMode(1);
for(x=0; x<nResults; x++) {
	X = getResult("XM", x);
		Y = getResult("YM", x);
				selectWindow("FDCM1 421.tif"); 
				Intensity= getPixel(X, Y); 
				selectWindow("WP heatmap centre.tif"); 
				CentreDistance= getPixel(X, Y); 
				selectWindow("WP DisEdge heat map.tif"); 
				EdgeDistance= getPixel(X, Y); 
				RatioCtoE= CentreDistance/EdgeDistance;
				print(Intensity,",", CentreDistance, ",", EdgeDistance,",", RatioCtoE); 
}//for x
setBatchMode(0);

waitForUser("Save values", "3. Save white pulp intenssity, distance centre/edge + c:e ratio values for graphpad analysis");

//***MZ***

//MZ WP dis FDC map

imageCalculator("AND create 32-bit", "FDC mask.tif","MZ heatmap.tif"); 
selectWindow("Result of FDC mask.tif"); 
run("35 step"); 
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/MZ WP FDC centre heatmap.tif"); 

imageCalculator("AND create 32-bit", "FDC mask.tif","MZ DisEdge heat map.tif"); 
selectWindow("Result of FDC mask.tif"); 
run("35 step"); 
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/MZ WP FDC edge heatmap.tif"); 

//MZ map dis from MZ into WP

selectWindow("MOMA mask.tif");
//newImage("DisEdge", "8-bit black", 1024, 1024, 1);
//roiManager("Select", MZ);
//run("Clear Outside");
//run("Make Inverse"); 
//changeValues(0, 254, 255);
run("Select None");
run("Distance Map");
run("35 step");
roiManager("Select", outline);
run("Clear Outside");
run("Make Inverse"); 
changeValues(0, 65535, 0);
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/WP dis from MZ heat map.tif");

print("\\Clear");
getRawStatistics(nPixels, mean, min, max, std, histogramArray); // random distribution control
for (n = 0; n < histogramArray.length; n++) {
	print(n,",",histogramArray [n]);
}
print(nPixels);

waitForUser("Save values", "4. Save white pulp edge random distribution control values for graphpad analysis");


newImage("DisEdge", "8-bit black", 1024, 1024, 1);
roiManager("Select", MZ);
changeValues(0, 254, 255);
run("Distance Map");
run("35 step");
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/MZ edge heat map.tif");

print("\\Clear");
getRawStatistics(nPixels, mean, min, max, std, histogramArray); // random distribution control
for (n = 0; n < histogramArray.length; n++) {
	print(n,",",histogramArray [n]);
}
print(nPixels);

waitForUser("Save values", "5. Save MZ centre random distribution control values for graphpad analysis");

//MZ dis FDC map
imageCalculator("AND create 32-bit", "FDC mask.tif","MZ DisEdge heat map.tif"); 
selectWindow("Result of FDC mask.tif"); 
run("35 step"); 
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/MZ FDC edge heatmap.tif"); 

//measure histo pixel distances

open("C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/FDCM1 421.tif");

selectWindow("FDC WP DisEdge heat map.tif");
getRawStatistics(nPixels, mean, min, max, std, dhistogram); 
selectWindow("FDC mask.tif");
selectWindow("Results")
print("nResults", nResults);
print("\\Clear");

setBatchMode(1);
for(x=0; x<nResults; x++) {
	X = getResult("XM", x);
		Y = getResult("YM", x);
				selectWindow("FDCM1 421.tif"); 
				IntensityWP= getPixel(X, Y); 
				selectWindow("MZ WP FDC centre heatmap.tif"); 
				CentreDistanceWP= getPixel(X, Y); 
				selectWindow("MZ WP FDC edge heatmap.tif"); 
				MZwpEdgeDistanceWP= getPixel(X, Y); 
				selectWindow("WP dis from MZ heat map.tif"); 
				EdgeDistanceWP= getPixel(X, Y); 
				selectWindow("MZ edge heat map.tif"); 
				MZEdgeDistanceWP= getPixel(X, Y); 
				print(IntensityWP,",", CentreDistanceWP, ",", EdgeDistanceWP,",", MZwpEdgeDistanceWP,",", MZEdgeDistanceWP); 
}//for x
setBatchMode(0);

waitForUser("Save values", "6. Save white pulp intenssity, distance WP centre/edge , distance from MZ into WP, c:e ratio values, distance MZ centre/edge for graphpad analysis");

//***B CELL FOLLICLES***

//Identrify white pulp region (convex hull)

 selectImage(ImageIDoriginal);
selectWindow("B220 mask.tif");
run("8-bit"); 
run("Options...", "iterations=10 count=1 edm=8-bit do=Open");
//run("Fill Holes"); 
run("Options...", "iterations=3 count=1 black edm=8-bit do=Dilate"); 
roiManager("Select", line);
run("Make Inverse");
changeValues(0, 65535, 0);
run("Analyze Particles...", "size=1000-Infinity display summarize add");

setMinAndMax(0, 255);
run("Select None");
run("Analyze Particles...", "size=100-Infinity show=[Count Masks] display summarize add");
run("glasbey");
run("Enhance Contrast", "saturated=0.35");

setOption("BlackBackground", true);
run("Convert to Mask");
run("Options...", "iterations=10 count=1 black edm=8-bit do=[Fill Holes]");

//run("Point Tool...");
//setTool("multipoint");
//waitForUser("Select B cell follicle")
//getSelectionCoordinates(xpoints, ypoints)

///for (o=0; o<xpoints.length; o=o+1){
	//value=getPixel(xpoints[o], ypoints[o]);
	//changeValues(value, value, 255);
	//print(o, value);
//}

//changeValues(0, 254, 0);
//changeValues(256,65000,0);
//run("Select None");
run("8-bit");
run("Options...", "iterations=1 count=1 black edm=8-bit do=Nothing");
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/Follicle mask.tif");
run("Duplicate...", "FolDuplicate");
rename("FolDuplicate");
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/Follicle mask-1.tif");

selectWindow("Follicle mask.tif");
run("8-bit");
run("Select None");
run("Grays");
setAutoThreshold("Huang dark");
run("Convert to Mask");
//run("Invert");
run("Create Selection");
roiManager("Add");
Follicle= roiManager("count")-1;
roiManager("Select", Follicle);
roiManager("Rename", "Follicle"); 
run("Make Binary");
run("Make Inverse");
run("Convex Hull");
run("Create Mask");
run("Create Selection");
roiManager("Add");
lineF= roiManager("count")-1;
roiManager("Select", lineF);
roiManager("Rename", "Follicle outline"); 
run("Clear Outside");
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/Follicle mask.tif");

//Identify centre

Roi.getContainedPoints(xpoints, ypoints);
Array.getStatistics(xpoints, min, max, meanx, stdDev);
Array.getStatistics(ypoints, min, max, meany, stdDev);
print(meanx, meany);

//Create EDM for heat map (centre)

getDimensions(width, height, channel, slices, frames); 
newImage("conditional dilation", "8-bit black", width, height, 1);
CondiFlID=getImageID();
setPixel(meanx, meany, 255); 
run("Options...", "iterations=3 count=1 edm=16-bit do=Nothing");
run("Distance Map");
run("35 step");
roiManager("Select", Follicle);
run("Make Inverse");
changeValues(0, 65535, 0);
run("Make Inverse");
getRawStatistics(nPixels, mean, min, max, std, histogram); 
setMinAndMax(min, max);
print(min,max);
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/Follicle heatmap centre.tif") 
print("\\Clear");
getRawStatistics(nPixels, mean, min, max, std, histogramArray); // random distribution control
for (n = 0; n < histogramArray.length; n++) {
	print(n,",",histogramArray [n]);
}
print(nPixels);

waitForUser("Save values", "7. Save B cell Follicle centre random distribution control values for graphpad analysis");

//FDC centre map

selectWindow("FDC mask.tif");
imageCalculator("AND create 32-bit", "FDC mask.tif","Follicle heatmap centre.tif"); 
selectWindow("Result of FDC mask.tif"); 
run("35 step"); 
roiManager("Select", Follicle);
run("Clear Outside");
run("Make Inverse");  
changeValues(0, 255, 0);
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/Follicle FDC heatmap.tif"); 

//Ki67 centre map

selectWindow("Ki67 mask.tif");
imageCalculator("AND create 32-bit", "Ki67 mask.tif","Follicle heatmap centre.tif"); 
selectWindow("Result of Ki67 mask.tif"); 
run("35 step"); 
roiManager("Select", Follicle);
run("Clear Outside");
run("Make Inverse");  
changeValues(0, 255, 0);
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/Follicle Ki67 heatmap.tif"); 

//Create EDM for heat map (distance)

newImage("DisEdge", "8-bit black", 1024, 1024, 1);
roiManager("Select", Follicle);
run("Clear Outside");
run("Make Inverse"); 
changeValues(0, 254, 255);
run("Distance Map");
run("35 step");
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/Follicle disEdge heat map.tif");

print("\\Clear");
getRawStatistics(nPixels, mean, min, max, std, histogramArray); // random distribution control
for (n = 0; n < histogramArray.length; n++) {
	print(n,",",histogramArray [n]);
}
print(nPixels);

waitForUser("Save values", "8. Save B cell follicle edge random distribution control values for graphpad analysis");

//FDC edge map

selectWindow("FDC mask.tif");
roiManager("Select", Follicle);
run("Clear Outside");
run("Make Inverse");
changeValues(0, 255, 0);
imageCalculator("AND create 32-bit", "FDC mask.tif","Follicle disEdge heat map.tif"); 
selectWindow("Result of FDC mask.tif"); 
run("35 step"); 
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/Follicle FDC disEdge heat map.tif"); 

//Ki67 edge map

selectWindow("Ki67 mask.tif");
imageCalculator("AND create 32-bit", "Ki67 mask.tif","Follicle disEdge heat map.tif"); 
selectWindow("Result of Ki67 mask.tif"); 
run("35 step"); 
roiManager("Select", Follicle);
run("Clear Outside");
run("Make Inverse");  
changeValues(0, 255, 0);
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/Follicle Ki67 disEdge heatmap.tif"); 

//measure histo pixel distances

open("C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/FDCM1 421.tif");

selectWindow("Follicle FDC disEdge heat map.tif");
getRawStatistics(nPixels, mean, min, max, std, dhistogram); 
selectWindow("FDC mask.tif");
selectWindow("Results")
print("nResults", nResults);
print("\\Clear");

setBatchMode(1);
for(x=0; x<nResults; x++) {
	X = getResult("XM", x);
		Y = getResult("YM", x);
				selectWindow("FDCM1 421.tif"); 
				IntensityF= getPixel(X, Y); 
				selectWindow("Follicle heatmap centre.tif"); 
				CentreDistanceF= getPixel(X, Y); 
				selectWindow("Follicle disEdge heat map.tif"); 
				EdgeDistanceF= getPixel(X, Y); 
				//print(IntensityF,",", CentreDistanceF, ",", EdgeDistanceF); 
}//for x
setBatchMode(0);

waitForUser("Save values", "9. Save B cell follicle FDC intenssity, distance centre/edge + c:e ratio values for graphpad analysis");

//***CLUSTER ANALYSIS***

nData= 8; //sets c loop for random controls below ie number of random controls
newImage("Data", "32-bit black", nData+1, 255, 1);

selectWindow("FDC mask.tif");
run("Select None");
run("Distance Map");
run("35 step");
roiManager("Select", Follicle);
run("Make Inverse");
changeValues(0, 65535, 0);
run("Make Inverse");
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/Follicle FDC Cluster heat map.tif"); 
getRawStatistics(nPixels, mean, min, max, std, GChistArray);
run("Select None");

selectWindow("Data");
areaGC=0;// cumulate area
areaRand=0;
for (v = 0; v < 255; v++) { // based on data image size
	if (v< GChistArray.length) areaGC= areaGC+GChistArray[v];
	print(areaGC);
	setPixel(0, v, areaGC);
} // creates image with data for biological data = always column 0, controls 1 etc

//RandomJ create loop with number of FDC ROIs as loop variable ie create same number of random dots in GC ROI

selectWindow("FDC mask.tif");
Height= getHeight();
Width= getWidth();

newImage("Random", "8-bit black", Height, Width, 1);
run("Add...", "value=1");
roiManager("Select", Follicle);
changeValues(0, 255, 0);
run("Select None");
RandomID=getImageID();

//RandomJID=getImageID();

selectWindow("FDC mask.tif");
preFDCno= roiManager("count");
run("Select None");
setThreshold(1, 255);
run("Analyze Particles...", "size=10-Infinity display add");
postFDCno= roiManager("count");
FDCno=postFDCno-preFDCno;
print("FDCno=", FDCno);
print("\\Clear");

for (c = 0; c < nData; c++) { 
noWanted=FDCno;
noHave=0;

selectImage(RandomID);
run("Duplicate...", "RandomDuplicate");
rename("RandomDuplicate");
Random=getImageID();

newImage("RandomJ", "16-bit black", 256, 2, 1);
run("RandomJ Uniform", "min=0.0 max=1023 insertion=Additive"); 
RandomJID=getImageID();

setBatchMode(1);
for (i = 0; i < 1024; i++) {
	selectImage(RandomJID);
	xpos=getValue(i, 0);
	ypos=getValue(i, 1);
	selectWindow("FDC mask.tif");
	roiManager("Select", Follicle);
	posLoc= selectionContains(xpos, ypos);
	if(posLoc==1) {
	print(i,xpos,ypos,posLoc);	
	selectImage(Random);
	roiWant= preFDCno+noHave+1;
	print("roiWant", roiWant);
	roiManager("select", roiWant); 
	getSelectionBounds(x, y, ROIwidth, ROIheight);
	Xpos= xpos-ROIwidth/2;
	Ypos= ypos+ROIheight/2;
	Roi.move(Xpos,Ypos);
	getRawStatistics(nPixels, mean, min, ROImax, std, histogram);
	if (ROImax==0)	fill();
	noHave=noHave+1;
	} // end if  i loop
	if(noHave==noWanted-1) break;
	}// end break loop

selectImage(Random);
run("Select None");
changeValues(1, 1, 0);
run("Distance Map");
run("35 step");
roiManager("Select", Follicle);
run("Make Inverse");
changeValues(0, 65535, 0);
SaveName= "Random cluster" + c;
SaveLoc= "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/" + SaveName + ".tif";
saveAs("Tiff", SaveLoc); 
roiManager("Select", Follicle);

RhistArray = "RhistArray"+c;
getRawStatistics(nPixels, mean, min, max, std, RhistArray);
print(c," array size ",RhistArray.length);
run("Select None");

iMax= maxOf(GChistArray.length, RhistArray.length);

selectWindow("Data");
areaGC=0;// cumulate area
areaRand=0;
for (v = 0; v < 255; v++) { //255 based on Data image size 
	//areaGC= areaGC+RhistArray[v];
	if (v< RhistArray.length) areaGC= areaGC+RhistArray[v];
	//print(v," val", val);
setPixel(c+1, v, areaGC);
} // creates image with data for all c loops

} // end of c loop

selectWindow("Data");
saveAs("Text Image", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/10 Data- Follicle cumulative area.txt");
setBatchMode(0);

//***GERMINAL CENTRES***

//Identify centre and create EDM for heat map

selectWindow("Ki67 mask.tif");
roiManager("Select", Follicle);
run("Make Inverse");
changeValues(0, 65535, 0);
run("Select None");
run("Duplicate...", " ");
roiManager("Select", Follicle);
run("Clear Outside");
run("Make Inverse");  
changeValues(0, 255, 0);
run("Select None");
run("Options...", "iterations=4 count=5 black edm=16-bit do=Erode");
GCIDoriginal=getImageID();
run("Duplicate...", "duplicate"); 
GCIDduplicate=getImageID();

run("Options...", "iterations=10 count=5 black edm=8-bit do=[Fill Holes]");
run("Options...", "iterations=10 count=11 black edm=8-bit do=[Erode]");
run("Options...", "iterations=10 count=11 black edm=8-bit do=[Dilate]");

setTool("polygon");
waitForUser("Outline GC ki67+ GC regions");
roiManager("Add");
GCOutlineSel= roiManager("count")-1;
roiManager("Select", GCOutlineSel);
roiManager("Rename", "GCOutlineSel");
setForegroundColor(255, 255, 255);
run("Fill", "slice");
roiManager("Select", GCOutlineSel);
run("Make Inverse");
changeValues(0, 65535, 0);
run("Select None");
run("Create Selection");
roiManager("Add");
GCOutline= roiManager("count")-1;
roiManager("Select", GCOutline);
roiManager("Rename", "GCOutline");
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/Ki67 GC mask.tif") 

run("Options...", "iterations=10 count=1 black edm=8-bit do=Erode");
preGCno= roiManager("count");
run("Analyze Particles...", "size=2000-Infinity show=[Count Masks] display summarize add");
PostGCcountImage= getImageID();
postGCno= roiManager("count");
GCno=postGCno-preGCno;
print(GCno, postGCno, preGCno);

getDimensions(width, height, channel, slices, frames); 
newImage("GC conditional dilation", "8-bit black", width, height, 1);

noGCWanted=GCno;
noGCHave=0;
print(noGCHave, noGCWanted); 

for (l = 0; l < noGCWanted; l++) { 

print("noGCHave=", noGCHave);

	selectImage(PostGCcountImage);
	roiGCWant= preGCno+noGCHave; 
	print("roiGCWant", roiGCWant);
	roiManager("select", roiGCWant);
	
	Roi.getContainedPoints(xpoints, ypoints);
	Array.getStatistics(xpoints, min, max, meanGCx, stdDev);
	Array.getStatistics(ypoints, min, max, meanGCy, stdDev);
	print(meanGCx, meanGCy);

	selectWindow("GC conditional dilation");
	setPixel(meanGCx, meanGCy, 255);

	noGCHave= noGCHave+1;

	print(noGCHave);
	
} // end of l loop

run("Options...", "iterations=3 count=1 edm=16-bit do=Nothing");
run("Distance Map");
run("35 step");
roiManager("Select", Follicle);
run("Make Inverse");
changeValues(0, 65535, 0);
run("Make Inverse");
getRawStatistics(nPixels, mean, min, max, std, histogram); 
setMinAndMax(min, max);
print(min,max);
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/Follicle GC centre heatmap.tif") 

print("\\Clear");
getRawStatistics(nPixels, mean, min, max, std, histogramArray); // random distribution control
for (n = 0; n < histogramArray.length; n++) {
	print(n,",",histogramArray [n]);
}
print(nPixels);

waitForUser("Save values", "11. Save Follicle GC centre random distribution control values for graphpad analysis");

run("Select None");
run("Duplicate...", "duplicate");
roiManager("Select", GCOutline);
run("Make Inverse");
changeValues(0, 65535, 0);
run("Make Inverse");
getRawStatistics(nPixels, mean, min, max, std, histogram); 
setMinAndMax(min, max);
print(min,max);
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/GC centre heatmap.tif") 

print("\\Clear");
getRawStatistics(nPixels, mean, min, max, std, histogramArray); // random distribution control
for (n = 0; n < histogramArray.length; n++) {
	print(n,",",histogramArray [n]);
}
print(nPixels);

waitForUser("Save values", "12. Save GC centre random distribution control values for graphpad analysis");


//FDC centre map

selectWindow("FDC mask.tif");
imageCalculator("AND create 32-bit", "FDC mask.tif","Follicle GC centre heatmap.tif"); 
selectWindow("Result of FDC mask.tif"); 
run("35 step"); 
roiManager("Select", Follicle);
run("Clear Outside");
run("Make Inverse");  
changeValues(0, 255, 0);
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/Follicle GC centre FDC heatmap.tif"); 

selectWindow("FDC mask.tif");
imageCalculator("AND create 32-bit", "FDC mask.tif","GC centre heatmap.tif"); 
selectWindow("Result of FDC mask.tif"); 
run("35 step"); 
roiManager("Select", Follicle);
run("Clear Outside");
run("Make Inverse");  
changeValues(0, 255, 0);
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/GC FDC heatmap.tif"); 

//Inside/outside GC distance maps

selectImage(GCIDduplicate);
run("Select None");
run("Duplicate...", "duplicate");
GCIDoutside=getImageID();
newImage("GCinside", "8-bit black", 1024, 1024, 1);
GCIDinside=getImageID();

selectImage(GCIDoutside);
run("Select None");
run("Distance Map");
run("35 step");
roiManager("Select", Follicle);
run("Make Inverse");
changeValues(0, 65535, 0);
run("Make Inverse");
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/Follicle GC outside heat map.tif"); 
print("\\Clear");
getRawStatistics(nPixels, mean, min, max, std, histogramArray); // random distribution control
for (n = 0; n < histogramArray.length; n++) {
	print(n,",",histogramArray [n]);
}
print(nPixels);

waitForUser("Save values", "13.Save GC-Follicle outside random distribution control values for graphpad analysis");


selectImage(GCIDinside);
roiManager("Select", GCOutline);
run("Make Inverse"); 
changeValues(0, 254, 255);
run("Distance Map");
run("35 step");
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/Follicle GC inside heat map.tif");

print("\\Clear");
getRawStatistics(nPixels, mean, min, max, std, histogramArray); // random distribution control
for (n = 0; n < histogramArray.length; n++) {
	print(n,",",histogramArray [n]);
}
print(nPixels);

waitForUser("Save values", "14 Save GC-Follicle inside random distribution control values for graphpad analysis");

//FDC centre maps

selectWindow("FDC mask.tif");
imageCalculator("AND create 32-bit", "FDC mask.tif","Follicle GC outside heat map.tif"); 
selectWindow("Result of FDC mask.tif"); 
run("35 step"); 
roiManager("Select", Follicle);
run("Clear Outside");
run("Make Inverse");  
changeValues(0, 255, 0);
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/Follicle GC outside FDC heat map.tif"); 

selectWindow("FDC mask.tif");
imageCalculator("AND create 32-bit", "FDC mask.tif","Follicle GC inside heat map.tif"); 
selectWindow("Result of FDC mask.tif"); 
run("35 step"); 
roiManager("Select", Follicle);
run("Clear Outside");
run("Make Inverse");  
changeValues(0, 255, 0);
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/Follicle GC inside FDC heat map.tif"); 

//measure histo pixel distances

selectWindow("FDC mask.tif");
getRawStatistics(nPixels, mean, min, max, std, dhistogram); 
selectWindow("FDC mask.tif");
selectWindow("Results")
print("nResults", nResults);
print("\\Clear");

setBatchMode(1);
for(x=0; x<nResults; x++) {
	X = getResult("XM", x);
		Y = getResult("YM", x);
				selectWindow("FDCM1 421.tif"); 
				IntensityGC= getPixel(X, Y); 
				selectWindow("Follicle GC centre heatmap.tif"); 
				GCCentreDistanceGC= getPixel(X, Y);
				selectWindow("GC centre heatmap.tif"); 
				CentreDistanceGC= getPixel(X, Y); 
				selectWindow("Follicle GC inside FDC heat map.tif");
				InsideDistanceGC= getPixel(X, Y); 
				selectWindow("Follicle GC outside FDC heat map.tif");
				OutsideDistanceGC= getPixel(X, Y); 
				print(IntensityGC,",", CentreDistanceGC,",", GCCentreDistanceGC,",", InsideDistanceGC, ",", OutsideDistanceGC); 
}//for x
setBatchMode(0);

waitForUser("Save values", "15. Save GC intenssity, fol/GC distance centre + c:e ratio values, inside/outside distance GC values for graphpad analysis");

//***CLUSTER ANALYSIS***

nData= 8; //sets c loop for random controls below ie number of random controls
newImage("Data", "32-bit black", nData+1, 255, 1);

selectWindow("FDC mask.tif");
run("Select None");
run("Distance Map");
run("35 step");
roiManager("Select", GCOutline);
run("Make Inverse");
changeValues(0, 65535, 0);
run("Make Inverse");
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/GC FDC Cluster heat map.tif"); 
getRawStatistics(nPixels, mean, min, max, std, GChistArray);
run("Select None");

selectWindow("Data");
areaGC=0;// cumulate area
areaRand=0;
for (v = 0; v < 255; v++) { // based on data image size
	if (v< GChistArray.length) areaGC= areaGC+GChistArray[v];
	print(areaGC);
	setPixel(0, v, areaGC);
} // creates image with data for biological data = always column 0, controls 1 etc

//RandomJ create loop with number of FDC ROIs as loop variable ie create same number of random dots in GC ROI

selectWindow("FDC mask.tif");
Height= getHeight();
Width= getWidth();

newImage("Random", "8-bit black", Height, Width, 1);
run("Add...", "value=1");
roiManager("Select", GCOutline);
changeValues(0, 255, 0);
run("Select None");
RandomID=getImageID();

//RandomJID=getImageID();

selectWindow("FDC mask.tif");
preFDCno= roiManager("count");
run("Select None");
setThreshold(1, 255);
run("Analyze Particles...", "size=10-Infinity display add");
postFDCno= roiManager("count");
FDCno=postFDCno-preFDCno;
print("FDCno=", FDCno);
print("\\Clear");

for (c = 0; c < nData; c++) { 
noWanted=FDCno;
noHave=0;

selectImage(RandomID);
run("Duplicate...", "RandomDuplicate");
rename("RandomDuplicate");
Random=getImageID();

newImage("RandomJ", "16-bit black", 256, 2, 1);
run("RandomJ Uniform", "min=0.0 max=1023 insertion=Additive"); 
RandomJID=getImageID();

setBatchMode(1);
for (i = 0; i < 1024; i++) {
	selectImage(RandomJID);
	xpos=getValue(i, 0);
	ypos=getValue(i, 1);
	selectWindow("FDC mask.tif");
	roiManager("Select", GCOutline);
	posLoc= selectionContains(xpos, ypos);
	if(posLoc==1) {
	print(i,xpos,ypos,posLoc);	
	selectImage(Random);
	roiWant= preFDCno+noHave+1;
	print("roiWant", roiWant);
	roiManager("select", roiWant); 
	getSelectionBounds(x, y, ROIwidth, ROIheight);
	Xpos= xpos-ROIwidth/2;
	Ypos= ypos+ROIheight/2;
	Roi.move(Xpos,Ypos);
	getRawStatistics(nPixels, mean, min, ROImax, std, histogram);
	if (ROImax==0)	fill();
	noHave=noHave+1;
	} // end if  i loop
	if(noHave==noWanted-1) break;
	}// end break loop

selectImage(Random);
run("Select None");
changeValues(1, 1, 0);
run("Distance Map");
run("35 step");
roiManager("Select", GCOutline);
run("Make Inverse");
changeValues(0, 65535, 0);
SaveName= "GC Random cluster" + c;
SaveLoc= "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/" + SaveName + ".tif";
saveAs("Tiff", SaveLoc); 
roiManager("Select", GCOutline);

RhistArray = "RhistArray"+c;
getRawStatistics(nPixels, mean, min, max, std, RhistArray);
print(c," array size ",RhistArray.length);
run("Select None");

iMax= maxOf(GChistArray.length, RhistArray.length);

selectWindow("Data");
areaGC=0;// cumulate area
areaRand=0;
for (v = 0; v < 255; v++) { //255 based on Data image size 
	//areaGC= areaGC+RhistArray[v];
	if (v< RhistArray.length) areaGC= areaGC+RhistArray[v];
	//print(v," val", val);
setPixel(c+1, v, areaGC);
} // creates image with data for all c loops

} // end of c loop

selectWindow("Data");
saveAs("Text Image", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/16 Data- GC cumulative area.txt");
setBatchMode(0);

//***LZ***
//Identify centre and create EDM for heat map

selectWindow("Ki67 mask.tif");
roiManager("Select", GCOutline);
run("Make Inverse");
changeValues(0, 65535, 0);
run("Duplicate...", "duplicate"); 
LZIDtemplate=getImageID();

setTool("polygon");
waitForUser("Outline GC ki67+ LZ regions");
roiManager("Add");
LZOutlineSel= roiManager("count")-1;
roiManager("Select", LZOutlineSel);
roiManager("Rename", "LZOutlineSel");
run("Make Inverse");
changeValues(0, 65535, 0);
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/LZ mask.tif");

run("Duplicate...", "LZ mask duplicate");
roiManager("Select", LZOutlineSel);
setForegroundColor(255, 255, 255);
run("Fill", "slice");
roiManager("Select", LZOutlineSel);
run("Make Inverse");
changeValues(0, 65535, 0);
run("Select None");
run("Invert");
run("Create Selection");
roiManager("Add");
LZOutline= roiManager("count")-1;
roiManager("Select", LZOutline);
roiManager("Rename", "LZOutline");

//Identify centre

Roi.getContainedPoints(xpoints, ypoints);
Array.getStatistics(xpoints, min, max, meanLZx, stdDev);
Array.getStatistics(ypoints, min, max, meanLZy, stdDev);
print(meanLZx, meanLZy);

//Create EDM for heat map (centre)

selectWindow("LZ mask-1.tif");
getDimensions(width, height, channel, slices, frames); 
newImage("LZconditional dilation", "8-bit black", width, height, 1);
LZcondilID=getImageID();
setPixel(meanLZx, meanLZy, 255); 

run("Options...", "iterations=3 count=1 edm=16-bit do=Nothing");
run("Distance Map");
run("35 step");
roiManager("Select", Follicle);
run("Make Inverse");
changeValues(0, 65535, 0);
run("Make Inverse");
getRawStatistics(nPixels, mean, min, max, std, histogram); 
setMinAndMax(min, max);
print(min,max);
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/Follicle LZ centre heatmap.tif") 

print("\\Clear");
getRawStatistics(nPixels, mean, min, max, std, histogramArray); // random distribution control
for (n = 0; n < histogramArray.length; n++) {
	print(n,",",histogramArray [n]);
}
print(nPixels);

waitForUser("Save values", "17. Save Follicle LZ centre random distribution control values for graphpad analysis");

run("Select None");
run("Duplicate...", "duplicate"); 
roiManager("Select", GCOutline);
run("Make Inverse");
changeValues(0, 65535, 0);
run("Make Inverse");
getRawStatistics(nPixels, mean, min, max, std, histogram); 
setMinAndMax(min, max);
print(min,max);
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/GC LZ centre heatmap.tif") 

print("\\Clear");
getRawStatistics(nPixels, mean, min, max, std, histogramArray); // random distribution control
for (n = 0; n < histogramArray.length; n++) {
	print(n,",",histogramArray [n]);
}
print(nPixels);

waitForUser("Save values", "18. Save GC LZ centre random distribution control values for graphpad analysis");


run("Select None");
run("Duplicate...", "duplicate");
roiManager("Select", LZOutline);
run("Make Inverse");
changeValues(0, 65535, 0);
run("Make Inverse");
getRawStatistics(nPixels, mean, min, max, std, histogram); 
setMinAndMax(min, max);
print(min,max);
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/LZ centre heatmap.tif") 

print("\\Clear");
getRawStatistics(nPixels, mean, min, max, std, histogramArray); // random distribution control
for (n = 0; n < histogramArray.length; n++) {
	print(n,",",histogramArray [n]);
}
print(nPixels);

waitForUser("Save values", "19. Save LZ centre random distribution control values for graphpad analysis");


//FDC centre map

selectWindow("FDC mask.tif");
imageCalculator("AND create 32-bit", "FDC mask.tif","Follicle LZ centre heatmap.tif");
selectWindow("Result of FDC mask.tif"); 
run("35 step"); 
roiManager("Select", Follicle);
run("Clear Outside");
run("Make Inverse");  
changeValues(0, 255, 0);
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/Follicle LZ centre FDC heatmap.tif"); 

run("35 step"); 
roiManager("Select", GCOutline);
run("Clear Outside");
run("Make Inverse");  
changeValues(0, 255, 0);
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/GC LZ centre FDC heatmap.tif"); 

roiManager("Select", LZOutline);
run("Clear Outside");
run("Make Inverse");  
changeValues(0, 255, 0);
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/LZ FDC heatmap.tif"); 

//Inside/outside LZ distance maps

newImage("LZoutside", "8-bit black", 1024, 1024, 1);
LZIDoutside=getImageID();
newImage("LZinside", "8-bit black", 1024, 1024, 1);
LZIDinside=getImageID();

selectImage(LZIDoutside);
roiManager("Select", LZOutline);
changeValues(0, 255, 255);
run("Select None");
run("Distance Map");
run("35 step");
roiManager("Select", Follicle);
run("Make Inverse");
changeValues(0, 65535, 0);
run("Make Inverse");
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/Follicle LZ outside heat map.tif"); 
print("\\Clear");
getRawStatistics(nPixels, mean, min, max, std, histogramArray); // random distribution control
for (n = 0; n < histogramArray.length; n++) {
	print(n,",",histogramArray [n]);
}
print(nPixels);

waitForUser("Save values", "20.Save LZ-Follicle outside random distribution control values for graphpad analysis");

run("Select None");
run("Duplicate...", "duplicate"); 
roiManager("Select", GCOutline);
run("Make Inverse");
changeValues(0, 65535, 0);
run("Make Inverse");
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/GC LZ outside heat map.tif"); 
print("\\Clear");
getRawStatistics(nPixels, mean, min, max, std, histogramArray); // random distribution control
for (n = 0; n < histogramArray.length; n++) {
	print(n,",",histogramArray [n]);
}
print(nPixels);

waitForUser("Save values", "21.Save LZ-GC outside random distribution control values for graphpad analysis");

selectImage(LZIDinside);
roiManager("Select", LZOutline);
run("Make Inverse"); 
changeValues(0, 254, 255);
run("Distance Map");
run("35 step");
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/Follicle LZ inside heat map.tif");

print("\\Clear");
getRawStatistics(nPixels, mean, min, max, std, histogramArray); // random distribution control
for (n = 0; n < histogramArray.length; n++) {
	print(n,",",histogramArray [n]);
}
print(nPixels);

waitForUser("Save values", "22 Save LZ-Follicle inside random distribution control values for graphpad analysis");

//FDC centre maps

selectWindow("FDC mask.tif");
imageCalculator("AND create 32-bit", "FDC mask.tif","Follicle LZ outside heat map.tif"); 
selectWindow("Result of FDC mask.tif"); 
run("35 step"); 
roiManager("Select", Follicle);
run("Clear Outside");
run("Make Inverse");  
changeValues(0, 255, 0);
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/Follicle LZ outside FDC heat map.tif"); 

selectWindow("FDC mask.tif");
imageCalculator("AND create 32-bit", "FDC mask.tif","Follicle LZ inside heat map.tif"); 
selectWindow("Result of FDC mask.tif"); 
run("35 step"); 
roiManager("Select", Follicle);
run("Clear Outside");
run("Make Inverse");  
changeValues(0, 255, 0);
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/Follicle LZ inside FDC heat map.tif"); 

//measure histo pixel distances

selectWindow("FDC mask.tif");
getRawStatistics(nPixels, mean, min, max, std, dhistogram); 
selectWindow("FDC mask.tif");
selectWindow("Results")
print("nResults", nResults);
print("\\Clear");

setBatchMode(1);
for(x=0; x<nResults; x++) {
	X = getResult("XM", x);
		Y = getResult("YM", x);
				selectWindow("FDCM1 421.tif"); 
				IntensityLZ= getPixel(X, Y); 
				selectWindow("Follicle LZ centre heatmap.tif"); 
				LZCentreDistanceLZ= getPixel(X, Y);
				selectWindow("LZ centre heatmap.tif"); 
				CentreDistanceLZ= getPixel(X, Y); 
				selectWindow("Follicle LZ inside FDC heat map.tif");
				InsideDistanceLZ= getPixel(X, Y); 
				selectWindow("Follicle LZ outside FDC heat map.tif");
				OutsideDistanceLZ= getPixel(X, Y); 
				print(IntensityLZ,",", CentreDistanceLZ,",", LZCentreDistanceLZ,",", InsideDistanceLZ, ",", OutsideDistanceLZ); 
}//for x
setBatchMode(0);

waitForUser("Save values", "23. Save LZ intenssity, fol/GC distance centre + c:e ratio values, inside/outside distance LZ values for graphpad analysis");

//***CLUSTER ANALYSIS***

nData= 8; //sets c loop for random controls below ie number of random controls
newImage("Data", "32-bit black", nData+1, 255, 1);

selectWindow("FDC mask.tif");
run("Select None");
run("Distance Map");
run("35 step");
roiManager("Select", LZOutline);
run("Make Inverse");
changeValues(0, 65535, 0);
run("Make Inverse");
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/LZ FDC Cluster heat map.tif"); 
getRawStatistics(nPixels, mean, min, max, std, LZhistArray);
run("Select None");

selectWindow("Data");
areaLZ=0;// cumulate area
areaRand=0;
for (v = 0; v < 255; v++) { // based on data image size
	if (v< LZhistArray.length) areaLZ= areaLZ+LZhistArray[v];
	print(areaLZ);
	setPixel(0, v, areaLZ);
} // creates image with data for biological data = always column 0, controls 1 etc

//RandomJ create loop with number of FDC ROIs as loop variable ie create same number of random dots in LZ ROI

selectWindow("FDC mask.tif");
Height= getHeight();
Width= getWidth();

newImage("Random", "8-bit black", Height, Width, 1);
run("Add...", "value=1");
roiManager("Select", LZOutline);
changeValues(0, 255, 0);
run("Select None");
RandomID=getImageID();

//RandomJID=getImageID();

selectWindow("FDC mask.tif");
preFDCno= roiManager("count");
run("Select None");
setThreshold(1, 255);
run("Analyze Particles...", "size=10-Infinity display add");
postFDCno= roiManager("count");
FDCno=postFDCno-preFDCno;
print("FDCno=", FDCno);
print("\\Clear");

for (c = 0; c < nData; c++) { 
noWanted=FDCno;
noHave=0;

selectImage(RandomID);
run("Duplicate...", "RandomDuplicate");
rename("RandomDuplicate");
Random=getImageID();

newImage("RandomJ", "16-bit black", 256, 2, 1);
run("RandomJ Uniform", "min=0.0 max=1023 insertion=Additive"); 
RandomJID=getImageID();

setBatchMode(1);
for (i = 0; i < 1024; i++) {
	selectImage(RandomJID);
	xpos=getValue(i, 0);
	ypos=getValue(i, 1);
	selectWindow("FDC mask.tif");
	roiManager("Select", LZOutline);
	posLoc= selectionContains(xpos, ypos);
	if(posLoc==1) {
	print(i,xpos,ypos,posLoc);	
	selectImage(Random);
	roiWant= preFDCno+noHave+1;
	print("roiWant", roiWant);
	roiManager("select", roiWant); 
	getSelectionBounds(x, y, ROIwidth, ROIheight);
	Xpos= xpos-ROIwidth/2;
	Ypos= ypos+ROIheight/2;
	Roi.move(Xpos,Ypos);
	getRawStatistics(nPixels, mean, min, ROImax, std, histogram);
	if (ROImax==0)	fill();
	noHave=noHave+1;
	} // end if  i loop
	if(noHave==noWanted-1) break;
	}// end break loop

selectImage(Random);
run("Select None");
changeValues(1, 1, 0);
run("Distance Map");
run("35 step");
roiManager("Select", LZOutline);
run("Make Inverse");
changeValues(0, 65535, 0);
SaveName= "LZ Random cluster" + c;
SaveLoc= "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/" + SaveName + ".tif";
saveAs("Tiff", SaveLoc); 
roiManager("Select", LZOutline);

RhistArray = "RhistArray"+c;
getRawStatistics(nPixels, mean, min, max, std, RhistArray);
print(c," array size ",RhistArray.length);
run("Select None");

iMax= maxOf(LZhistArray.length, RhistArray.length);

selectWindow("Data");
areaLZ=0;// cumulate area
areaRand=0;
for (v = 0; v < 255; v++) { //255 based on Data image size 
	//areaLZ= areaLZ+RhistArray[v];
	if (v< RhistArray.length) areaLZ= areaLZ+RhistArray[v];
	//print(v," val", val);
setPixel(c+1, v, areaLZ);
} // creates image with data for all c loops

} // end of c loop

selectWindow("Data");
saveAs("Text Image", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/24 Data- LZ cumulative area.txt");
setBatchMode(0);

//***DZ***
//Identify centre and create EDM for heat map

newImage("DZ mask", "8-bit black", 1024, 1024, 1);
DZmask=getImageID();
roiManager("Select", GCOutline);
changeValues(0, 255, 255);
roiManager("Select", LZOutline);
changeValues(0, 255, 0);
run("Options...", "iterations=10 count=1 edm=16-bit do=Erode");
run("Options...", "iterations=10 count=1 edm=16-bit do=Dilate");
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/DZ mask.tif");

run("Select None");
run("Invert");
run("Create Selection");
roiManager("Add");
DZOutline= roiManager("count")-1;
roiManager("Select", DZOutline);
roiManager("Rename", "DZOutline");

//Identify centre

Roi.getContainedPoints(xpoints, ypoints);
Array.getStatistics(xpoints, min, max, meanDZx, stdDev);
Array.getStatistics(ypoints, min, max, meanDZy, stdDev);
print(meanDZx, meanDZy);

//Create EDM for heat map (centre)

selectWindow("DZ mask.tif");
getDimensions(width, height, channel, slices, frames); 
newImage("DZconditional dilation", "8-bit black", width, height, 1);
DZcondilID=getImageID();
setPixel(meanDZx, meanDZy, 255); 

run("Options...", "iterations=3 count=1 edm=16-bit do=Nothing");
run("Distance Map");
run("35 step");
roiManager("Select", Follicle);
run("Make Inverse");
changeValues(0, 65535, 0);
run("Make Inverse");
getRawStatistics(nPixels, mean, min, max, std, histogram); 
setMinAndMax(min, max);
print(min,max);
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/Follicle DZ centre heatmap.tif") 

print("\\Clear");
getRawStatistics(nPixels, mean, min, max, std, histogramArray); // random distribution control
for (n = 0; n < histogramArray.length; n++) {
	print(n,",",histogramArray [n]);
}
print(nPixels);

waitForUser("Save values", "25. Save Follicle DZ centre random distribution control values for graphpad analysis");

run("Select None");
run("Duplicate...", "duplicate"); 
roiManager("Select", GCOutline);
run("Make Inverse");
changeValues(0, 65535, 0);
run("Make Inverse");
getRawStatistics(nPixels, mean, min, max, std, histogram); 
setMinAndMax(min, max);
print(min,max);
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/GC DZ centre heatmap.tif") 

print("\\Clear");
getRawStatistics(nPixels, mean, min, max, std, histogramArray); // random distribution control
for (n = 0; n < histogramArray.length; n++) {
	print(n,",",histogramArray [n]);
}
print(nPixels);

waitForUser("Save values", "26. Save GC DZ centre random distribution control values for graphpad analysis");


run("Select None");
run("Duplicate...", "duplicate");
roiManager("Select", DZOutline);
run("Make Inverse");
changeValues(0, 65535, 0);
run("Make Inverse");
getRawStatistics(nPixels, mean, min, max, std, histogram); 
setMinAndMax(min, max);
print(min,max);
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/DZ centre heatmap.tif") 

print("\\Clear");
getRawStatistics(nPixels, mean, min, max, std, histogramArray); // random distribution control
for (n = 0; n < histogramArray.length; n++) {
	print(n,",",histogramArray [n]);
}
print(nPixels);

waitForUser("Save values", "27. Save DZ centre random distribution control values for graphpad analysis");


//FDC centre map

selectWindow("FDC mask.tif");
imageCalculator("AND create 32-bit", "FDC mask.tif","Follicle DZ centre heatmap.tif");
selectWindow("Result of FDC mask.tif"); 
run("35 step"); 
roiManager("Select", Follicle);
run("Clear Outside");
run("Make Inverse");  
changeValues(0, 255, 0);
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/Follicle DZ centre FDC heatmap.tif"); 

run("35 step"); 
roiManager("Select", GCOutline);
run("Clear Outside");
run("Make Inverse");  
changeValues(0, 255, 0);
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/GC DZ centre FDC heatmap.tif"); 

roiManager("Select", DZOutline);
run("Clear Outside");
run("Make Inverse");  
changeValues(0, 255, 0);
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/DZ FDC heatmap.tif"); 

//Inside/outside DZ distance maps

newImage("GCinside", "8-bit black", 1024, 1024, 1);
DZIDoutside=getImageID();
newImage("GCinside", "8-bit black", 1024, 1024, 1);
DZIDinside=getImageID();

selectImage(DZIDoutside);
roiManager("Select", DZOutline);
changeValues(0, 255, 255);
run("Select None");
run("Distance Map");
run("35 step");
roiManager("Select", Follicle);
run("Make Inverse");
changeValues(0, 65535, 0);
run("Make Inverse");
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/Follicle DZ outside heat map.tif"); 
print("\\Clear");
getRawStatistics(nPixels, mean, min, max, std, histogramArray); // random distribution control
for (n = 0; n < histogramArray.length; n++) {
	print(n,",",histogramArray [n]);
}
print(nPixels);

waitForUser("Save values", "28.Save DZ-Follicle outside random distribution control values for graphpad analysis");

run("Select None");
run("Duplicate...", "duplicate"); 
roiManager("Select", GCOutline);
run("Make Inverse");
changeValues(0, 65535, 0);
run("Make Inverse");
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/GC DZ outside heat map.tif"); 
print("\\Clear");
getRawStatistics(nPixels, mean, min, max, std, histogramArray); // random distribution control
for (n = 0; n < histogramArray.length; n++) {
	print(n,",",histogramArray [n]);
}
print(nPixels);

waitForUser("Save values", "28.Save DZ-GC outside random distribution control values for graphpad analysis");

selectImage(DZIDinside);
roiManager("Select", DZOutline);
run("Make Inverse"); 
changeValues(0, 254, 255);
run("Distance Map");
run("35 step");
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/Follicle DZ inside heat map.tif");

print("\\Clear");
getRawStatistics(nPixels, mean, min, max, std, histogramArray); // random distribution control
for (n = 0; n < histogramArray.length; n++) {
	print(n,",",histogramArray [n]);
}
print(nPixels);

waitForUser("Save values", "29 Save DZ-Follicle inside random distribution control values for graphpad analysis");

//FDC centre maps

selectWindow("FDC mask.tif");
imageCalculator("AND create 32-bit", "FDC mask.tif","Follicle DZ outside heat map.tif"); 
selectWindow("Result of FDC mask.tif"); 
run("35 step"); 
roiManager("Select", Follicle);
run("Clear Outside");
run("Make Inverse");  
changeValues(0, 255, 0);
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/Follicle DZ outside FDC heat map.tif"); 

selectWindow("FDC mask.tif");
imageCalculator("AND create 32-bit", "FDC mask.tif","Follicle DZ inside heat map.tif"); 
selectWindow("Result of FDC mask.tif"); 
run("35 step"); 
roiManager("Select", Follicle);
run("Clear Outside");
run("Make Inverse");  
changeValues(0, 255, 0);
saveAs("Tiff", "C:/Users/jessa/Desktop/GC Zone Macro/Macro Output/Follicle DZ inside FDC heat map.tif"); 

//measure histo pixel distances

selectWindow("FDC mask.tif");
getRawStatistics(nPixels, mean, min, max, std, dhistogram); 
selectWindow("FDC mask.tif");
selectWindow("Results")
print("nResults", nResults);
print("\\Clear");

setBatchMode(1);
for(x=0; x<nResults; x++) {
	X = getResult("XM", x);
		Y = getResult("YM", x);
				selectWindow("FDCM1 421.tif"); 
				IntensityDZ= getPixel(X, Y); 
				selectWindow("Follicle DZ centre heatmap.tif"); 
				DZCentreDistancDZe= getPixel(X, Y);
				selectWindow("DZ centre heatmap.tif"); 
				CentreDistanceDZ= getPixel(X, Y); 
				selectWindow("Follicle DZ inside FDC heat map.tif");
				InsideDistanceDZ= getPixel(X, Y); 
				selectWindow("Follicle DZ outside FDC heat map.tif");
				OutsideDistanceDZ= getPixel(X, Y); 
				print(IntensityDZ,",", CentreDistanceDZ,",", DZCentreDistanceDZ,",", InsideDistanceDZ, ",", OutsideDistanceDZ); 
}//for x
setBatchMode(0);

waitForUser("Save values", "30. Save DZ intenssity, fol/GC distance centre + c:e ratio values, inside/outside distance DZ values for graphpad analysis");

print("\\Clear");


run("Close");
run("Close All"); 
print("\\Close"); 
roiManager("Reset"); 
