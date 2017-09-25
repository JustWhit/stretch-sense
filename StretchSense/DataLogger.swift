import Foundation

class DataLogger: NSObject {
    
    var counter:Int = 0;
    var file: FileHandle!;
    
    override init() {
        super.init();
        do {
            let fileNameUrl = URL(string: FileManager.default.currentDirectoryPath + "/stretchsense.log");
            file = try FileHandle(forWritingTo: fileNameUrl!);
            file.write("pF,time,sample\n".data(using: String.Encoding.utf8, allowLossyConversion: false)!);
        } catch {
            print(error.localizedDescription);
            exit(EXIT_FAILURE);
        }
    }
    
    func writeStretchSenseData(stretchSenseEntry : String) {
        file.seekToEndOfFile();
        let toWrite = stretchSenseEntry + "\n";
        file.write(toWrite.data(using: String.Encoding.utf8, allowLossyConversion: false)!);
    }
    
    deinit {
        file.closeFile();
    }
    
}

