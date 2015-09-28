import src.FileIO.SiemensReader as SiemensReader
import src.MRI.util.processing_utils as processing_utils


sr = SiemensReader.SiemensReader(path)
f = sr.read()
