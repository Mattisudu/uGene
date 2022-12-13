class OrdMap {
    constructor() {
        this.keys = [];
        this.values = [];
    }
    // Method
    set(key, value) {
        if(this.keys.slice(-1)[0] < key){
            this.keys.push(key);
            this.values.push(value)
            return true
        }
        if(this.keys[0]>key){
            this.keys.unshift(key);
            this.values.unshift(value)
            return true

        }

        let b_l = 0
        let b_r = this.keys.length
        let p = Math.floor(this.keys.length/2)
        while(b_l!==b_r) {
            if (this.keys[p] === key) return false;
            if (key > this.keys[p]){b_l = p + 1; p = b_l + Math.floor((b_r - b_l) / 2)}
            else{b_r = p; p = b_l + Math.floor((b_r - b_l) / 2);}
        }
        this.keys.splice(p,0,key);
        this.values.splice(p,0,value);
        return true;
    }
    setDebug(key, value) {
        this.keys.push(key)
        this.keys.push(value)
        return true;
    }
    setOverride(key, value) {
        if(this.keys.slice(-1)[0] < key){
            this.keys.push(key);
            this.values.push(value)
            return true
        }
        if(this.keys[0]>key){
            this.keys.unshift(key);
            this.values.unshift(value)
            return true

        }

        let b_l = 0
        let b_r = this.keys.length
        let p = Math.floor(this.keys.length/2)
        while(b_l!==b_r) {
            if (this.keys[p] === key){
                this.keys.splice(p,1,key);
                this.values.splice(p,1,value);
                return true;
            }
            if (key > this.keys[p]){b_l = p + 1; p = b_l + Math.floor((b_r - b_l) / 2)}
            else{b_r = p; p = b_l + Math.floor((b_r - b_l) / 2);}
        }
        this.keys.splice(p,0,key);
        this.values.splice(p,0,value);
        return true;
    }
    size(){return this.keys.length;}
    [Symbol.iterator]() {
        var index = -1;
        return {next: () => ({ value: this.values[++index], done: !(index in this.values) })};
    }
}
