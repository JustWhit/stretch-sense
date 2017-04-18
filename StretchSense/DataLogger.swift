import Foundation

class DataLogger: NSObject {
    
    var counter:Int = 0;
    
    override init() {
        super.init();
    }
    
    func writeToFile(data: String) {
        do {
            let fileName = FileManager.default.currentDirectoryPath + "/stretchsense.log";
            let fileNameUrl = URL(string: fileName);
            let file = try FileHandle(forWritingTo: fileNameUrl!);
            file.seekToEndOfFile();
            let date = Date();
            let formatter = DateFormatter();
            formatter.dateFormat = "HH:mm:ss";
            let ts = formatter.string(from: date);
            let toWrite = data + "," + ts + ",\(counter)" + "\n";
            file.write(toWrite.data(using: String.Encoding.utf8, allowLossyConversion: false)!);
            file.closeFile();
            counter += 1;
        } catch {
            print(error.localizedDescription);
            exit(EXIT_FAILURE);
        }
    }
    
}

