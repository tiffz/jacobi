// Basic MatrixWrapper class
function MatrixWrapper(r, c) {
  if ((r instanceof Array) && (r[0] instanceof Array)) {
    this.matrix = r;
    this.rows = r.length;
    this.cols = r[0].length;
  } else {
    this.rows = r;
    this.cols = c;
    if (!this.rows) {
      // Default value
      this.rows = 1;
    }
    if (!this.cols) {
      this.cols = this.rows;
    }

    this.matrix = new Array(this.rows);
    for (var i = 0; i < this.rows; i++) {
      this.matrix[i] = new Array(this.cols);
    }
  }

  this.set = function(r, c, val) {
    this.matrix[r][c] = val;
  }

  this.get = function(r, c) {
    if (this.matrix[r][c] == null || isNaN(this.matrix[r][c])) {
      this.matrix[r][c] = 0;
    }
    return this.matrix[r][c];
  }

  this.getArray = function() {
    var array = new Array(this.rows);
    for (var i = 0; i < this.rows; i++) {
      array[i] = new Array(this.cols);
      for (var j = 0; j < this.cols; j++) {
        array[i][j] = this.get(i, j);
      }
    }
    return array;
  }

  this.s = function() {
    this.sylvester = $M(this.getArray());
    return this.sylvester;
  }

  this.toString = function(delimiter) {
    if (!delimiter) {
      delimiter = ' ';
    }
    var str = '';
    for (var i = 0; i < this.rows; i++) {
      for (var j = 0; j < this.cols; j++) {
        str = str + this.get(i, j) + delimiter;
      }
      str = str + '\n';
    }
    return str;
  }

  this.getOffB = function() {
    var offB = 0;
    for (var i = 0; i < this.rows; i++) {
      for (var j = i + 1; j < this.cols; j++) {
        offB += Math.pow(this.get(i, j), 2);
        offB += Math.pow(this.get(j, i), 2);
      }
    }
    return offB;
  }

  // One iteration of Jacobi's algorithm
  // Matrix MUST be symmetric!
  this.iterate = function(r, c) {
    var str = '';
    if (!isNumber(r) || !isNumber(c)) {
      var offD = this.getLargestOffDiagonal();
      r = offD[0];
      c = offD[1];
    }
    var sub = this.getSubmatrix(r, c);

    var a = sub.get(0, 0);
    var b = sub.get(0, 1);
    var d = sub.get(1, 1);

    var next = null;
    if (b == 0) {
      next = this;
    } else {

      // Get the largest eigenvalue
      var e1 = (a + d) / 2 + Math.sqrt(Math.pow((a - d) / 2, 2) + (b * b));

      // B is a Sylvester
      var B = sub.subtract(e1);


      // Get the first row of the matrix
      var r1 = B.row(1);


      // Get orthogonal vector
      // Note Sylvester indices start at 1
      var u1 = $V([-r1.e(2), r1.e(1)]).toUnitVector();
      var u2 = $V([u1.e(2), -u1.e(1)]);

      // Get the final U
      var u = $M([[u1.e(1), u2.e(1)],
                  [u1.e(2), u2.e(2)]]
                );
      var iN = Matrix.I(this.rows);
      var G = iN.elements;
      G[r][r] = u1.e(1);
      G[r][c] = u1.e(2);
      G[c][r] = u2.e(1);
      G[c][c] = u2.e(2);
      G = $M(G);

      var GT = G.transpose();
      var Gn = GT.multiply(this.s()).multiply(G);
      var next = new MatrixWrapper(Gn.elements);
    }
    return next;
  }

  this.subtract = function(value) {
    if (value instanceof MatrixWrapper) {
      value = value.s();
    } else if (isNumber(value)) {
      value = Matrix.I(this.rows).multiply(value);
    }
    return this.s().subtract(value);
  }

  this.getLargestOffDiagonal = function() {
    var r = 0;
    var c = 1;
    var max = Math.abs(this.get(r, c));
    for (var i = 0; i < this.rows; i++) {
      for (var j = i + 1; j < this.cols; j++) {
        var num = Math.abs(this.get(i, j));
        if (num > max) {
          max = num;
          r = i;
          c = j;
        }
      }
    }

    return [r, c, max];
  }

  this.getSubmatrix = function(r, c) {
    var sub = [[this.get(r, r), this.get(r, c)], [this.get(c, r), this.get(c, c)]];
    sub = new MatrixWrapper(sub);
    return sub;
  }

  this.prettyPrint = function() {
    var matrixClass = 'matrix';
    var base = 50;
    var bracket_size = 50;
    if (this.rows > 9) {
      base = 20;
      matrixClass = 'big_matrix';
      bracket_size = 80;
    }
    var height = base * this.rows;
    var width = base * this.cols + bracket_size * 2;
    bracket_style = 'style="height: ' + height + 'px; line-height: ' + height + 'px; font-size:' + (height + 30) + 'px;"';
    str = '<div class="' + matrixClass + '" style="height: ' + height + 'px; width: ' + width + 'px;">';
    str = str + '<div class="bracket left" ' + bracket_style + '>[</div>';
    str = str + '<div class="bracket right" ' + bracket_style + '>[</div>';
    for (var i = 0; i < this.rows; i++) {
      for (var j = 0; j < this.cols; j++) {
        var num = Math.round(this.get(i, j) * 1000) / 1000;
        if (Math.abs(num) > 1000) {
          num = num.toExponential();
        }
        str = str + '<div class="matrix_entry">' + num + '</div>';
      }
      str = str + '\n';
    }
    str += '</div>';
    return str;
  }

  this.makeRandomSymmetric = function(min, max) {
    if (this.rows != this.cols) {
      return;
    }
    if (!min || !max) {
      min = -10;
      max = 10;
    }
    for (var i = 0; i < this.rows; i++) {
      for (var j = i; j < this.cols; j++) {
        var num = Math.floor(Math.random() * (max - min + 1)) + min;
        this.set(i, j, num);
        this.set(j, i, num);
      }
    }
  }

  this.isSymmetric = function() {
    var symmetric = true;
    if (this.rows != this.cols) {
      symmetric = false;
    }

    var i = 0;
    while (symmetric == true && i < this.rows) {
      var j = i + 1;
      while (symmetric == true && j < this.cols) {
        if (this.get(i, j) != this.get(j, i)) {
          symmetric = false;
        }
        j++;
      }
      i++;
    }
    return symmetric
  }
}

function textToMatrixWrapper(text, delim) {
  if (!delim) {
    delim = ' ';
  }
  var lines = text.split("\n");
  var m = [];
  var k = 0;
  for (var i = 0; i < lines.length; i++) {
    lines[i] = lines[i].split(delim);
    var l = 0;
    for (var j = 0; j < lines[i].length; j++) {
      var num = parseFloat(lines[i][j]) * 1;
      if (num != undefined && !isNaN(num)) {
        if (l == 0) {
          // Initialize 2D array
          m[k] = [];
        }
        m[k][l] = num;
        l++;
      }
    }
    if (l > 0) {

      k++;
    }
  }

  return new MatrixWrapper(m);
}

function isNumber (o) {
  return !isNaN (o-0) && o !== null && o !== false;
}