//preparation des fichiers tiff

sVer = 1.22;
saveDir = getDir("choix du repertoire de sauvegarde");
list = getList("image.titles");
if (list.length==0)
	print("aucune image ouverte");
	
else {
    //sauvegarde et fermeture
    for (i=0; i<list.length; i++){
		print("  "+list[i]);
		selectWindow(list[i]);
		saveAs("tiff", saveDir + list[i] + ".tiff");
		close();
		open(saveDir + list[i] + ".tiff", "virtual");
		temp=getImageID();
		Stack.getDimensions(width, height, channels, slices, frames);
		if (frames>1){
			setBatchMode("hide");
			j=1;
			while (j<=frames){
				Stack.setFrame(j);
				getRawStatistics(nPixels, mean, min, max, std, histogram);
				if (mean==0 && max ==0 && min==0){
					run("Delete Slice");
					j=j-1;
				}
				j++;
				Stack.getDimensions(width, height, channels, slices, frames);
			}
			setBatchMode("exit and display");
			selectImage(temp);
			wait(2);
			if (bitDepth()==16) run("8-bit");
			wait(2);
			saveAs("tiff", saveDir + list[i] + ".tiff");
			wait(2);
			close();
		}
	}
}
