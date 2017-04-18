import Foundation

let runLoop = RunLoop.current;
let distantFuture = Date.distantFuture;

let worker = BluetoothWorker();

while (runLoop.run(mode: RunLoopMode.defaultRunLoopMode, before: distantFuture)){
    print("### listening for a stretchsense device");
}
