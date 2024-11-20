from ErrorDetection import *

if __name__ == "__main__":
    ioe = "Data/PRECAI.DAT"
    gcp = "Data/PRECKI.DAT"
    pixel = "Data/PREPHIy10.DAT"
    output = "Data/result.txt"

    test = ErrorDetection()
    test.process(ioe, gcp, pixel, output)
