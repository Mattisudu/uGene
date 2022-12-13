let time_start = undefined;
let last_stage = undefined;

function timerReset(restart = undefined){
    time_start = undefined;
    last_stage = undefined;

    if (restart !== undefined){
        startTimer()
    }
}
function startTimer(){
    time_start = Date.now()
    last_stage = time_start
}
function lastTimerReset(){
    if(time_start === undefined) {
        startTimer()
    }
    else{
        last_stage = Date.now()
    }
}
function logRuntime(message=""){
    if(time_start === undefined){
        console.log("First run startTimer()")
        startTimer()
    }
    else{
    const now = Date.now()
    console.log("Passing time = ", now - last_stage, " Total passing time = ", now - time_start," ", message);
    last_stage = now

    }
}

function performTest2(func, args){
    const  strart = Date.now()
    func(args);
    const end = Date.now()
    console.log("performTest() result: ", end-strart,"    ", func)
    return end - strart
}