TODO add GPU edit distance support for sequence modified myers distance.

Since VRAM is becoming much larger in response to rise of AI needs, it may be faster to do distance calculation and bk-tree searches on GPU. In my use cases, the edit distance calculation is incredibly quick, its the moving the data to CPU is the main bottle neck. Doing this part on GPU would up CPU that I can use for follow up steps without impeding the bk-tree searches.
