// ---------- Helper functions ----------

// Convert a single character to digit
function charToDigit(ch) {
  const code = ch.toLowerCase().charCodeAt(0);
  if (ch >= '0' && ch <= '9') return code - 48;
  if (ch >= 'a' && ch <= 'z') return 10 + (code - 97);
  throw new Error("Invalid digit " + ch);
}

// Convert string from given base to BigInt
function convertFromBase(str, base) {
  let val = 0n;
  for (let ch of str) {
    let d = BigInt(charToDigit(ch));
    if (d >= BigInt(base)) throw new Error(`Digit '${ch}' >= base ${base}`);
    val = val * BigInt(base) + d;
  }
  return val;
}

// GCD for Fraction reduction
function gcd(a, b) {
  while (b !== 0n) {
    [a, b] = [b, a % b];
  }
  return a;
}

// Fraction class for exact arithmetic
class Fraction {
  constructor(n, d = 1n) {
    if (d === 0n) throw new Error("Divide by zero");
    if (d < 0n) { n = -n; d = -d; }
    const g = gcd(n < 0n ? -n : n, d);
    this.n = n / g;
    this.d = d / g;
  }
  add(o) { return new Fraction(this.n * o.d + o.n * this.d, this.d * o.d); }
  sub(o) { return new Fraction(this.n * o.d - o.n * this.d, this.d * o.d); }
  mul(o) { return new Fraction(this.n * o.n, this.d * o.d); }
  div(o) { return new Fraction(this.n * o.d, this.d * o.n); }
  equals(o) { return this.n === o.n && this.d === o.d; }
  toString() { return this.d === 1n ? this.n.toString() : this.n + "/" + this.d; }
}

// Lagrange interpolation
function lagrangeInterpolate(points) {
  const k = points.length;
  let coeffs = Array(k).fill(null).map(() => new Fraction(0n,1n));

  for (let j = 0; j < k; j++) {
    let [xj, yj] = points[j];
    let numerator = [new Fraction(1n,1n)];
    let denom = new Fraction(1n,1n);

    for (let m = 0; m < k; m++) {
      if (m === j) continue;
      let [xm] = points[m];
      numerator = numerator.map((c, idx) => {
        if (idx === 0) return c.mul(new Fraction(-xm,1n));
        return c.mul(new Fraction(-xm,1n)).add(numerator[idx-1]);
      });
      numerator.push(numerator[numerator.length-2] || new Fraction(1n,1n));
      denom = denom.mul(new Fraction(xj - xm, 1n));
    }

    for (let d = 0; d < numerator.length; d++) {
      coeffs[d] = coeffs[d].add(numerator[d].mul(new Fraction(yj,1n)).div(denom));
    }
  }
  return coeffs;
}

// Evaluate polynomial at x
function evalPoly(coeffs, x) {
  let sum = new Fraction(0n,1n);
  let xp = new Fraction(1n,1n);
  for (let c of coeffs) {
    sum = sum.add(c.mul(xp));
    xp = xp.mul(new Fraction(x,1n));
  }
  return sum;
}

// Generate k-combinations
function* kCombinations(arr, k) {
  const n = arr.length;
  function* helper(start, combo) {
    if (combo.length === k) {
      yield combo;
      return;
    }
    for (let i = start; i < n; i++) {
      yield* helper(i+1, combo.concat([arr[i]]));
    }
  }
  yield* helper(0, []);
}

// ---------- Main solver ----------
function solveShares(input){
  const {keys} = input;
  const k = keys.k;
  const shares = [];

  // If provided_shares exists, only take those
  const allowed = input.provided_shares || Object.keys(input).filter(k => k !== "keys" && k !== "provided_shares");

  for (let key of allowed) {
    if (key === "keys" || key === "provided_shares") continue;
    const x = BigInt(key);
    const base = parseInt(input[key].base, 10);
    const val = convertFromBase(input[key].value, base);
    shares.push([x, val]);
  }

  if (shares.length < k) {
    throw new Error(`Not enough shares provided (need at least ${k})`);
  }

  let best = null;
  for (const subset of kCombinations(shares, k)) {
    const coeffs = lagrangeInterpolate(subset);
    const inliers = [], outliers = [];
    for (const [x,y] of shares) {
      const yhat = evalPoly(coeffs, x);
      if (yhat.equals(new Fraction(y,1n))) inliers.push([x.toString(), y.toString()]);
      else outliers.push([x.toString(), y.toString()]);
    }
    if (!best || inliers.length > best.inliers.length) best = { coeffs, inliers, outliers };
  }
  if (!best) throw new Error("No polynomial found");

  return {
    secret: best.coeffs[0].toString(),
    wrong_shares: best.outliers,
    used_shares: best.inliers,
    polynomial: Object.fromEntries(best.coeffs.map((c,i)=>[`a${i}`,c.toString()])),
    notes: `Reconstructed degree-${k-1} polynomial`
  };
}

// ---------- Example runs ----------

// Small test
const input1 = {
  "keys": { "n": 4, "k": 3 },
  "1": { "base": "10", "value": "4" },
  "2": { "base": "2", "value": "111" },
  "3": { "base": "10", "value": "12" },
  "6": { "base": "4", "value": "213" },
  "provided_shares": ["1", "2", "3"]
};

// Big test with 10 shares
const input2 = {
  "keys": { "n": 10, "k": 7 },
  "1": { "base": "6", "value": "13444211440455345511" },
  "2": { "base": "15", "value": "aed7015a346d635" },
  "3": { "base": "15", "value": "6aeeb69631c227c" },
  "4": { "base": "16", "value": "e1b5e05623d881f" },
  "5": { "base": "8", "value": "316034514573652620673" },
  "6": { "base": "3", "value": "2122212201122002221120200210011020220200" },
  "7": { "base": "3", "value": "20120221122211000100210021102001201112121" },
  "8": { "base": "6", "value": "20220554335330240002224253" },
  "9": { "base": "12", "value": "45153788322a1255483" },
  "10": { "base": "7", "value": "1101613130313526312514143" }
};

// Choose which input to run
console.log("Small test:");
console.log(JSON.stringify(solveShares(input1), null, 2));

console.log("\nBig test:");
console.log(JSON.stringify(solveShares(input2), null, 2));
