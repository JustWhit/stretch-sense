import Foundation

class DataLogger: NSObject {
    
    var counter:Int = 0;
    var file: FileHandle!;
    let formatter = DateFormatter();
    
    override init() {
        super.init();
        formatter.dateFormat = "yyyy-MM-dd'T'HH:mm:ss.SSS";
        let date = String(format:"%f", (NSDate.timeIntervalSinceReferenceDate));
        do {
            // determine if the log file exists
            let fileManager = FileManager.default;
            let fileNameUrl = URL(string: fileManager.currentDirectoryPath + "/Data/CAP_" + date + ".csv");
            if (fileManager.fileExists(atPath: fileNameUrl!.path)) {
                do {
                    try fileManager.removeItem(atPath: fileNameUrl!.path)
                } catch {
                    print("ERROR:: unable to remove old capacitance data");
                    print(error.localizedDescription);
                    exit(EXIT_FAILURE);
                }
            }
            
            fileManager.createFile(atPath: fileNameUrl!.path, contents: nil, attributes: nil);
            file = try FileHandle(forWritingTo: fileNameUrl!);
            file.write("pF,time\n".data(using: String.Encoding.utf8, allowLossyConversion: false)!);
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

