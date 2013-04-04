// Generated by CoffeeScript 1.6.2
(function() {
  var dot, grad, perm,
    __indexOf = [].indexOf || function(item) { for (var i = 0, l = this.length; i < l; i++) { if (i in this && this[i] === item) return i; } return -1; },
    __slice = [].slice;

  grad = [[1.0, 1.0, 0.0], [-1.0, 1.0, 0.0], [1.0, -1.0, 0.0], [-1.0, -1.0, 0.0], [1.0, 0.0, 1.0], [-1.0, 0.0, 1.0], [1.0, 0.0, -1.0], [-1.0, 0.0, -1.0], [0.0, 1.0, 1.0], [0.0, -1.0, 1.0], [0.0, 1.0, -1.0], [0.0, -1.0, -1.0]];

  perm = [151, 160, 137, 91, 90, 15, 131, 13, 201, 95, 96, 53, 194, 233, 7, 225, 140, 36, 103, 30, 69, 142, 8, 99, 37, 240, 21, 10, 23, 190, 6, 148, 247, 120, 234, 75, 0, 26, 197, 62, 94, 252, 219, 203, 117, 35, 11, 32, 57, 177, 33, 88, 237, 149, 56, 87, 174, 20, 125, 136, 171, 168, 68, 175, 74, 165, 71, 134, 139, 48, 27, 166, 77, 146, 158, 231, 83, 111, 229, 122, 60, 211, 133, 230, 220, 105, 92, 41, 55, 46, 245, 40, 244, 102, 143, 54, 65, 25, 63, 161, 1, 216, 80, 73, 209, 76, 132, 187, 208, 89, 18, 169, 200, 196, 135, 130, 116, 188, 159, 86, 164, 100, 109, 198, 173, 186, 3, 64, 52, 217, 226, 250, 124, 123, 5, 202, 38, 147, 118, 126, 255, 82, 85, 212, 207, 206, 59, 227, 47, 16, 58, 17, 182, 189, 28, 42, 223, 183, 170, 213, 119, 248, 152, 2, 44, 154, 163, 70, 221, 153, 101, 155, 167, 43, 172, 9, 129, 22, 39, 253, 19, 98, 108, 110, 79, 113, 224, 232, 178, 185, 112, 104, 218, 246, 97, 228, 251, 34, 242, 193, 238, 210, 144, 12, 191, 179, 162, 241, 81, 51, 145, 235, 249, 14, 239, 107, 49, 192, 214, 31, 181, 199, 106, 157, 184, 84, 204, 176, 115, 121, 50, 45, 127, 4, 150, 254, 138, 236, 205, 93, 222, 114, 67, 29, 24, 72, 243, 141, 128, 195, 78, 66, 215, 61, 156, 180, 151, 160, 137, 91, 90, 15, 131, 13, 201, 95, 96, 53, 194, 233, 7, 225, 140, 36, 103, 30, 69, 142, 8, 99, 37, 240, 21, 10, 23, 190, 6, 148, 247, 120, 234, 75, 0, 26, 197, 62, 94, 252, 219, 203, 117, 35, 11, 32, 57, 177, 33, 88, 237, 149, 56, 87, 174, 20, 125, 136, 171, 168, 68, 175, 74, 165, 71, 134, 139, 48, 27, 166, 77, 146, 158, 231, 83, 111, 229, 122, 60, 211, 133, 230, 220, 105, 92, 41, 55, 46, 245, 40, 244, 102, 143, 54, 65, 25, 63, 161, 1, 216, 80, 73, 209, 76, 132, 187, 208, 89, 18, 169, 200, 196, 135, 130, 116, 188, 159, 86, 164, 100, 109, 198, 173, 186, 3, 64, 52, 217, 226, 250, 124, 123, 5, 202, 38, 147, 118, 126, 255, 82, 85, 212, 207, 206, 59, 227, 47, 16, 58, 17, 182, 189, 28, 42, 223, 183, 170, 213, 119, 248, 152, 2, 44, 154, 163, 70, 221, 153, 101, 155, 167, 43, 172, 9, 129, 22, 39, 253, 19, 98, 108, 110, 79, 113, 224, 232, 178, 185, 112, 104, 218, 246, 97, 228, 251, 34, 242, 193, 238, 210, 144, 12, 191, 179, 162, 241, 81, 51, 145, 235, 249, 14, 239, 107, 49, 192, 214, 31, 181, 199, 106, 157, 184, 84, 204, 176, 115, 121, 50, 45, 127, 4, 150, 254, 138, 236, 205, 93, 222, 114, 67, 29, 24, 72, 243, 141, 128, 195, 78, 66, 215, 61, 156, 180];

  dot = function(x, y, z, g) {
    return x * g[0] + y * g[1] + z * g[2];
  };

  window.noise = function(xin, yin, zin) {
    var F3, G3, X0, Y0, Z0, gi0, gi1, gi2, gi3, i, i1, i2, ii, j, j1, j2, jj, k, k1, k2, kk, n0, n1, n2, n3, s, t, t0, t1, t2, t3, x0, x1, x2, x3, y0, y1, y2, y3, z0, z1, z2, z3;

    F3 = 1.0 / 3.0;
    s = (xin + yin + zin) * F3;
    i = Math.floor(xin + s);
    j = Math.floor(yin + s);
    k = Math.floor(zin + s);
    G3 = 1.0 / 6.0;
    t = (i + j + k) * G3;
    X0 = i - t;
    Y0 = j - t;
    Z0 = k - t;
    x0 = xin - X0;
    y0 = yin - Y0;
    z0 = zin - Z0;
    if (x0 >= y0) {
      if (y0 >= z0) {
        i1 = 1;
        j1 = 0;
        k1 = 0;
        i2 = 1;
        j2 = 1;
        k2 = 0;
      } else if (x0 >= z0) {
        i1 = 1;
        j1 = 0;
        k1 = 0;
        i2 = 1;
        j2 = 0;
        k2 = 1;
      } else {
        i1 = 0;
        j1 = 0;
        k1 = 1;
        i2 = 1;
        j2 = 0;
        k2 = 1;
      }
    } else {
      if (y0 < z0) {
        i1 = 0;
        j1 = 0;
        k1 = 1;
        i2 = 0;
        j2 = 1;
        k2 = 1;
      } else if (x0 < z0) {
        i1 = 0;
        j1 = 1;
        k1 = 0;
        i2 = 0;
        j2 = 1;
        k2 = 1;
      } else {
        i1 = 0;
        j1 = 1;
        k1 = 0;
        i2 = 1;
        j2 = 1;
        k2 = 0;
      }
    }
    x1 = x0 - i1 + G3;
    y1 = y0 - j1 + G3;
    z1 = z0 - k1 + G3;
    x2 = x0 - i2 + 2.0 * G3;
    y2 = y0 - j2 + 2.0 * G3;
    z2 = z0 - k2 + 2.0 * G3;
    x3 = x0 - 1.0 + 3.0 * G3;
    y3 = y0 - 1.0 + 3.0 * G3;
    z3 = z0 - 1.0 + 3.0 * G3;
    ii = i & 255;
    jj = j & 255;
    kk = k & 255;
    gi0 = perm[ii + perm[jj + perm[kk]]] % 12;
    gi1 = perm[ii + i1 + perm[jj + j1 + perm[kk + k1]]] % 12;
    gi2 = perm[ii + i2 + perm[jj + j2 + perm[kk + k2]]] % 12;
    gi3 = perm[ii + 1 + perm[jj + 1 + perm[kk + 1]]] % 12;
    n0 = n1 = n2 = n3 = 0;
    t0 = 0.6 - x0 * x0 - y0 * y0 - z0 * z0;
    if (t0 > 0) {
      t0 *= t0;
      n0 = t0 * t0 * dot(x0, y0, z0, grad[gi0]);
    }
    t1 = 0.6 - x1 * x1 - y1 * y1 - z1 * z1;
    if (t1 > 0) {
      t1 *= t1;
      n1 = t1 * t1 * dot(x1, y1, z1, grad[gi1]);
    }
    t2 = 0.6 - x2 * x2 - y2 * y2 - z2 * z2;
    if (t2 > 0) {
      t2 *= t2;
      n2 = t2 * t2 * dot(x2, y2, z2, grad[gi2]);
    }
    t3 = 0.6 - x3 * x3 - y3 * y3 - z3 * z3;
    if (t3 > 0) {
      t3 *= t3;
      n3 = t3 * t3 * dot(x3, y3, z3, grad[gi3]);
    }
    return 16.0 * (n0 + n1 + n2 + n3) + 1.0;
  };

  window.simplex_noise = function(octaves, x, y, z) {
    var i, value, _i, _ref;

    value = 0.0;
    for (i = _i = 0, _ref = octaves - 1; 0 <= _ref ? _i <= _ref : _i >= _ref; i = 0 <= _ref ? ++_i : --_i) {
      value += noise(x * Math.pow(2, i), y * Math.pow(2, i), z * Math.pow(2, i));
    }
    return value;
  };

  $(function() {
    var calc_occs, camera, cube, cube_loop, doit, g, i, idx, idxa, idxt, idxvs, lookup, m, neigh, nsh, occs_per_node, occs_per_vert, place, precalc, r, rasterize, rays, renderer, s, scene, sphere_dist, to_canv, ts, _i, _len;

    renderer = new THREE.WebGLRenderer();
    camera = new THREE.PerspectiveCamera(45, 400 / 300, 0.1, 10000);
    scene = new THREE.Scene();
    renderer.setClearColor(new THREE.Color(0, 1));
    renderer.setSize(800, 600);
    $("body").append(renderer.domElement);
    s = 64;
    ts = 512;
    $('#canvas').attr({
      width: ts,
      height: ts
    });
    nsh = s / -2;
    idxt = function(x, y, z, ss) {
      ss || (ss = s);
      return x + y * ss + z * ss * ss;
    };
    idxa = [];
    (function() {
      var x, y, z, _i, _ref, _results;

      _results = [];
      for (x = _i = 0, _ref = s - 1; 0 <= _ref ? _i <= _ref : _i >= _ref; x = 0 <= _ref ? ++_i : --_i) {
        idxa[x] = [];
        _results.push((function() {
          var _j, _ref1, _results1;

          _results1 = [];
          for (y = _j = 0, _ref1 = s - 1; 0 <= _ref1 ? _j <= _ref1 : _j >= _ref1; y = 0 <= _ref1 ? ++_j : --_j) {
            idxa[x][y] = [];
            _results1.push((function() {
              var _k, _ref2, _results2;

              _results2 = [];
              for (z = _k = 0, _ref2 = s - 1; 0 <= _ref2 ? _k <= _ref2 : _k >= _ref2; z = 0 <= _ref2 ? ++_k : --_k) {
                _results2.push(idxa[x][y][z] = idxt(x, y, z));
              }
              return _results2;
            })());
          }
          return _results1;
        })());
      }
      return _results;
    })();
    idx = function(x, y, z) {
      return idxa[x][y][z];
    };
    cube_loop = function(fn) {
      var r, x, y, z, _i, _j, _k, _ref, _ref1, _ref2;

      r = [];
      for (z = _i = 0, _ref = s - 1; 0 <= _ref ? _i <= _ref : _i >= _ref; z = 0 <= _ref ? ++_i : --_i) {
        for (y = _j = 0, _ref1 = s - 1; 0 <= _ref1 ? _j <= _ref1 : _j >= _ref1; y = 0 <= _ref1 ? ++_j : --_j) {
          for (x = _k = 0, _ref2 = s - 1; 0 <= _ref2 ? _k <= _ref2 : _k >= _ref2; x = 0 <= _ref2 ? ++_k : --_k) {
            r.push(fn(x, y, z, idx(x, y, z)));
          }
        }
      }
      return r;
    };
    (function() {
      var v;

      v = new THREE.Vector3(nsh, nsh, nsh);
      return THREE.Vector3.prototype.repos = _.partial(THREE.Vector3.prototype.add, v);
    })();
    sphere_dist = function(n) {
      var dist, i, inc, o, phi, r, y, _i;

      dist = [];
      inc = Math.PI * (3 - Math.sqrt(5));
      o = 2 / n;
      for (i = _i = 1; 1 <= n ? _i <= n : _i >= n; i = 1 <= n ? ++_i : --_i) {
        y = i * o - 1 + (o / 2);
        r = Math.sqrt(1) - y * y;
        phi = i * inc;
        dist.push([Math.cos(phi) * r, y, Math.sin(phi) * r]);
      }
      return dist;
    };
    rasterize = function(x, y, z) {
      var cx, cy, cz, depth, ncx, ncy, ncz, points, sx, sy, sz;

      sx = x * 0.2;
      sy = y * 0.2;
      sz = z * 0.2;
      x = x < 0 ? -1 : 1;
      y = y < 0 ? -1 : 1;
      z = z < 0 ? -1 : 1;
      cx = 0;
      cy = 0;
      cz = 0;
      points = [];
      while (points.length < s) {
        ncx = Math.floor(x);
        ncy = Math.floor(y);
        ncz = Math.floor(z);
        if (ncx !== cx || ncy !== cy || ncz !== cz) {
          depth = Math.sqrt(x * x + y * y + z * z);
          cx = ncx;
          cy = ncy;
          cz = ncz;
          points.push(new THREE.Vector3(cx, cy, cz));
        }
        x += sx;
        y += sy;
        z += sz;
      }
      return points;
    };
    camera.position.z = s * 2.5;
    lookup = function(c, x, y, z) {
      if ((x >= 0 && x < s) && (y >= 0 && y < s) && (z >= 0 && z < s)) {
        return c[idx(x, y, z)];
      } else {
        return false;
      }
    };
    place = function(c, x, y, z, v) {
      var r, _i, _ref, _results;

      r = (function() {
        _results = [];
        for (var _i = 0, _ref = s - 1; 0 <= _ref ? _i <= _ref : _i >= _ref; 0 <= _ref ? _i++ : _i--){ _results.push(_i); }
        return _results;
      }).apply(this);
      if (__indexOf.call(r, x) >= 0 && __indexOf.call(r, y) >= 0 && __indexOf.call(r, z) >= 0) {
        return c[idx(x, y, z)] = v;
      } else {
        return false;
      }
    };
    neigh = _.memoize(function(cube, x, y, z) {
      var i, n, _i, _len, _results;

      n = [[x - 1, y, z], [x + 1, y, z], [x, y - 1, z], [x, y + 1, z], [x, y, z - 1], [x, y, z + 1]];
      _results = [];
      for (_i = 0, _len = n.length; _i < _len; _i++) {
        i = n[_i];
        _results.push(lookup.apply(null, [cube].concat(__slice.call(i))));
      }
      return _results;
    }, (function(c, x, y, z) {
      return [x, y, z].join();
    }));
    cube = cube_loop(function(x, y, z) {
      var caves, center_falloff, density, plateau_falloff, tv, tv3, tv5;

      tv = new THREE.Vector3(x, y, z).multiplyScalar(1 / s);
      tv3 = new THREE.Vector3(tv.x + 1, tv.y + 1, tv.z + 1).multiplyScalar(3);
      tv5 = tv.clone().multiplyScalar(5);
      if (tv.y <= 0.8) {
        plateau_falloff = 1.0;
      } else if (tv.y > 0.8 && tv.y < 0.9) {
        plateau_falloff = 1.0 - (tv.y - 0.8) * 10.0;
      } else {
        plateau_falloff = 0.0;
      }
      center_falloff = 0.1 / (Math.pow((tv.x - 0.5) * 1.5, 2) + Math.pow((tv.y - 1.0) * 0.8, 2) + Math.pow((tv.z - 0.5) * 1.5, 2));
      caves = Math.pow(simplex_noise(1, tv5.x, tv5.y, tv5.z), 3);
      density = simplex_noise(5, tv.x, tv.y / 2, tv.z) * plateau_falloff * center_falloff * Math.pow(simplex_noise(1, tv3.x, tv3.y, tv3.z) + 0.4, 1.8);
      if (caves < 0.5) {
        density = 0;
      }
      return density > 3.1;
    });
    if (s < 5) {
      console.log(cube);
    }
    rays = sphere_dist(50);
    for (_i = 0, _len = rays.length; _i < _len; _i++) {
      r = rays[_i];
      r.push(rasterize.apply(null, r), (function(func, args, ctor) {
        ctor.prototype = func.prototype;
        var child = new ctor, result = func.apply(child, args);
        return Object(result) === result ? result : child;
      })(THREE.Vector3, r, function(){}).normalize());
    }
    to_canv = function(inp) {
      var canvas, canvasHeight, canvasWidth, ctx, data, imageData, index, x, y, _j, _k, _ref, _ref1;

      canvas = document.getElementById('canvas');
      canvasWidth = canvas.width;
      canvasHeight = canvas.height;
      ctx = canvas.getContext('2d');
      imageData = ctx.getImageData(0, 0, canvasWidth, canvasHeight);
      data = imageData.data;
      for (y = _j = 0, _ref = canvasHeight - 1; 0 <= _ref ? _j <= _ref : _j >= _ref; y = 0 <= _ref ? ++_j : --_j) {
        for (x = _k = 0, _ref1 = canvasWidth - 1; 0 <= _ref1 ? _k <= _ref1 : _k >= _ref1; x = 0 <= _ref1 ? ++_k : --_k) {
          index = (y * canvasWidth + x) * 4;
          data[index] = inp[index];
          data[++index] = inp[index];
          data[++index] = inp[index];
          data[++index] = inp[index];
        }
      }
      return ctx.putImageData(imageData, 0, 0);
    };
    calc_occs = function() {
      var comp, cube_faces, derp, face, g, gl, i, input, l, l4, m, o, occ, output, ray, ray_index, rs, sc, sub, t, to_xyz, u, x, y, z, _j, _k, _l, _len1, _len2, _len3, _len4, _len5, _m, _n, _o, _ref;

      r = new THREE.WebGLRenderer({
        preserveDrawingBuffer: true
      });
      r.setSize(ts, ts);
      derp = [];
      for (_j = 0, _len1 = cube.length; _j < _len1; _j++) {
        sub = cube[_j];
        if (sub) {
          derp.push(255, 255, 255, 255);
        } else {
          derp.push(0, 0, 0, 255);
        }
      }
      for (i = _k = 1, _ref = ts * ts * 4 - derp.length; 1 <= _ref ? _k <= _ref : _k >= _ref; i = 1 <= _ref ? ++_k : --_k) {
        derp.push(0, 0, 0, 255);
      }
      input = new Uint8Array(derp);
      t = new THREE.DataTexture(input, ts, ts, THREE.RGBAFormat, THREE.UnsignedByteType, new THREE.UVMapping(), void 0, void 0, THREE.NearestFilter, THREE.NearestFilter);
      t.needsUpdate = true;
      rs = rays[0][3];
      u = {
        t: {
          type: 't',
          value: t
        },
        s: {
          type: 'f',
          value: s
        },
        ts: {
          type: 'f',
          value: ts
        },
        rayoffs: {
          type: 'v3v',
          value: rs
        }
      };
      sc = new THREE.Scene();
      g = new THREE.PlaneGeometry(2, 2);
      m = new THREE.Mesh(g, new THREE.ShaderMaterial({
        uniforms: u,
        vertexShader: $('#occv').text(),
        fragmentShader: $('#occf').text().replace('rayoffs[64]', "rayoffs[" + s + "]").replace('r < 64', "r < " + s)
      }));
      sc.add(m);
      output = new Uint8Array(ts * ts * 4);
      occ = new Array(ts * ts);
      to_xyz = function(cux, cuy) {
        var cox, coy, coz, fidx;

        fidx = cux + cuy * ts;
        coz = Math.floor(fidx / (s * s));
        fidx -= coz * s * s;
        coy = Math.floor(fidx / s);
        fidx -= coy * s;
        cox = fidx;
        return [cox, coy, coz];
      };
      gl = r.getContext();
      for (ray_index = _l = 0, _len2 = rays.length; _l < _len2; ray_index = ++_l) {
        ray = rays[ray_index];
        console.log(ray_index);
        u.rayoffs.value = ray[3];
        console.log("render gooo");
        r.render(sc, camera);
        console.log("render end");
        console.log("read gooo");
        gl.readPixels(0, 0, ts, ts, gl.RGBA, gl.UNSIGNED_BYTE, output);
        console.log("read end");
        console.log("total gooo");
        for (l = _m = 0, _len3 = output.length; _m < _len3; l = _m += 4) {
          comp = output[l];
          l4 = l / 4;
          if (ray_index === 0) {
            occ[l4] = [[0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0]];
          }
          x = ray[4].x;
          y = ray[4].y;
          z = ray[4].z;
          if (x < 0) {
            occ[l4][0][0] -= x;
            if (comp > 0) {
              occ[l4][0][1] -= x;
            }
          }
          if (x > 0) {
            occ[l4][1][0] += x;
            if (comp > 0) {
              occ[l4][1][1] += x;
            }
          }
          if (y < 0) {
            occ[l4][2][0] -= y;
            if (comp > 0) {
              occ[l4][2][1] -= y;
            }
          }
          if (y > 0) {
            occ[l4][3][0] += y;
            if (comp > 0) {
              occ[l4][3][1] += y;
            }
          }
          if (z < 0) {
            occ[l4][4][0] -= z;
            if (comp > 0) {
              occ[l4][4][1] -= z;
            }
          }
          if (z > 0) {
            occ[l4][5][0] += z;
            if (comp > 0) {
              occ[l4][5][1] += z;
            }
          }
        }
        console.log("total end");
      }
      o = [];
      for (_n = 0, _len4 = occ.length; _n < _len4; _n++) {
        cube_faces = occ[_n];
        t = [];
        for (_o = 0, _len5 = cube_faces.length; _o < _len5; _o++) {
          face = cube_faces[_o];
          t.push(face[0] && face[1] ? face[1] / face[0] : 0);
        }
        o.push(t);
      }
      return o;
    };
    occs_per_node = calc_occs();
    g = new THREE.Geometry();
    idxvs = function(x, y, z) {
      g.vertices.push(new THREE.Vector3(x, y, z).repos());
      return g.vertices.length - 1;
    };
    occs_per_vert = [];
    cube_loop(function(x, y, z, i) {
      if (cube[i]) {
        return _.each(neigh(cube, x, y, z), function(e2, i2) {
          var f, _, _j, _results;

          if (!e2[0]) {
            f = (function() {
              switch (i2) {
                case 0:
                  return new THREE.Face4(idxvs(x, y, z), idxvs(x, y, z + 1), idxvs(x, y + 1, z + 1), idxvs(x, y + 1, z));
                case 1:
                  return new THREE.Face4(idxvs(x + 1, y, z), idxvs(x + 1, y + 1, z), idxvs(x + 1, y + 1, z + 1), idxvs(x + 1, y, z + 1));
                case 2:
                  return new THREE.Face4(idxvs(x, y, z), idxvs(x + 1, y, z), idxvs(x + 1, y, z + 1), idxvs(x, y, z + 1));
                case 3:
                  return new THREE.Face4(idxvs(x, y + 1, z), idxvs(x, y + 1, z + 1), idxvs(x + 1, y + 1, z + 1), idxvs(x + 1, y + 1, z));
                case 4:
                  return new THREE.Face4(idxvs(x, y, z), idxvs(x, y + 1, z), idxvs(x + 1, y + 1, z), idxvs(x + 1, y, z));
                case 5:
                  return new THREE.Face4(idxvs(x, y, z + 1), idxvs(x + 1, y, z + 1), idxvs(x + 1, y + 1, z + 1), idxvs(x, y + 1, z + 1));
              }
            })();
            g.faces.push(f);
            _results = [];
            for (_ = _j = 1; _j <= 4; _ = ++_j) {
              _results.push(occs_per_vert.push(occs_per_node[idx(x, y, z)][i2]));
            }
            return _results;
          }
        });
      }
    });
    if (s < 5) {
      console.log(g);
    }
    m = new THREE.Mesh(g, new THREE.ShaderMaterial({
      vertexShader: $('#shlightv').text(),
      fragmentShader: $('#shlightf').text(),
      uniforms: {
        ScaleFactor: 1.0
      },
      attributes: {
        ocv: {
          type: 'f',
          value: occs_per_vert
        }
      }
    }));
    scene.add(m);
    g.computeFaceNormals();
    g.computeVertexNormals();
    g.normalsNeedUpdate = true;
    precalc = (function() {
      var _j, _ref, _results;

      _results = [];
      for (i = _j = 0.0, _ref = 2 * Math.PI; 0.0025 > 0 ? _j <= _ref : _j >= _ref; i = _j += 0.0025) {
        _results.push([s * Math.sin(i), s * Math.cos(i)]);
      }
      return _results;
    })();
    doit = function(t) {
      var u;

      u = precalc[Math.floor(t % precalc.length)] || [0, 1];
      camera.position.x = u[0];
      camera.position.y = u[1];
      camera.lookAt(scene.position);
      renderer.render(scene, camera);
      return requestAnimationFrame(doit);
    };
    return doit(0);
  });

}).call(this);
