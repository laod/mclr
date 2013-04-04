#Thanks to Florian Boesch http://codeflow.org/entries/2010/dec/09/minecraft-like-rendering-experiments-in-opengl-4/
#for the inspiration and significant chunks of code.

grad = [
  [1.0,1.0,0.0],[-1.0,1.0,0.0],[1.0,-1.0,0.0],[-1.0,-1.0,0.0],
  [1.0,0.0,1.0],[-1.0,0.0,1.0],[1.0,0.0,-1.0],[-1.0,0.0,-1.0],
  [0.0,1.0,1.0],[0.0,-1.0,1.0],[0.0,1.0,-1.0],[0.0,-1.0,-1.0]
]

perm = [151, 160, 137, 91, 90, 15, 131, 13, 201, 95, 96, 53, 194, 233, 7, 225, 140, 36, 103, 30, 69, 142, 8, 99, 37, 240, 21, 10, 23, 190, 6, 148, 247, 120, 234, 75, 0, 26, 197, 62, 94, 252, 219, 203, 117, 35, 11, 32, 57, 177, 33, 88, 237, 149, 56, 87, 174, 20, 125, 136, 171, 168, 68, 175, 74, 165, 71, 134, 139, 48, 27, 166, 77, 146, 158, 231, 83, 111, 229, 122, 60, 211, 133, 230, 220, 105, 92, 41, 55, 46, 245, 40, 244, 102, 143, 54, 65, 25, 63, 161, 1, 216, 80, 73, 209, 76, 132, 187, 208, 89, 18, 169, 200, 196, 135, 130, 116, 188, 159, 86, 164, 100, 109, 198, 173, 186, 3, 64, 52, 217, 226, 250, 124, 123, 5, 202, 38, 147, 118, 126, 255, 82, 85, 212, 207, 206, 59, 227, 47, 16, 58, 17, 182, 189, 28, 42, 223, 183, 170, 213, 119, 248, 152, 2, 44, 154, 163, 70, 221, 153, 101, 155, 167, 43, 172, 9, 129, 22, 39, 253, 19, 98, 108, 110, 79, 113, 224, 232, 178, 185, 112, 104, 218, 246, 97, 228, 251, 34, 242, 193, 238, 210, 144, 12, 191, 179, 162, 241, 81, 51, 145, 235, 249, 14, 239, 107, 49, 192, 214, 31, 181, 199, 106, 157, 184, 84, 204, 176, 115, 121, 50, 45, 127, 4, 150, 254, 138, 236, 205, 93, 222, 114, 67, 29, 24, 72, 243, 141, 128, 195, 78, 66, 215, 61, 156, 180, 151, 160, 137, 91, 90, 15, 131, 13, 201, 95, 96, 53, 194, 233, 7, 225, 140, 36, 103, 30, 69, 142, 8, 99, 37, 240, 21, 10, 23, 190, 6, 148, 247, 120, 234, 75, 0, 26, 197, 62, 94, 252, 219, 203, 117, 35, 11, 32, 57, 177, 33, 88, 237, 149, 56, 87, 174, 20, 125, 136, 171, 168, 68, 175, 74, 165, 71, 134, 139, 48, 27, 166, 77, 146, 158, 231, 83, 111, 229, 122, 60, 211, 133, 230, 220, 105, 92, 41, 55, 46, 245, 40, 244, 102, 143, 54, 65, 25, 63, 161, 1, 216, 80, 73, 209, 76, 132, 187, 208, 89, 18, 169, 200, 196, 135, 130, 116, 188, 159, 86, 164, 100, 109, 198, 173, 186, 3, 64, 52, 217, 226, 250, 124, 123, 5, 202, 38, 147, 118, 126, 255, 82, 85, 212, 207, 206, 59, 227, 47, 16, 58, 17, 182, 189, 28, 42, 223, 183, 170, 213, 119, 248, 152, 2, 44, 154, 163, 70, 221, 153, 101, 155, 167, 43, 172, 9, 129, 22, 39, 253, 19, 98, 108, 110, 79, 113, 224, 232, 178, 185, 112, 104, 218, 246, 97, 228, 251, 34, 242, 193, 238, 210, 144, 12, 191, 179, 162, 241, 81, 51, 145, 235, 249, 14, 239, 107, 49, 192, 214, 31, 181, 199, 106, 157, 184, 84, 204, 176, 115, 121, 50, 45, 127, 4, 150, 254, 138, 236, 205, 93, 222, 114, 67, 29, 24, 72, 243, 141, 128, 195, 78, 66, 215, 61, 156, 180]

dot = (x, y, z, g) ->
  x*g[0] + y*g[1] + z*g[2]

window.noise = (xin, yin, zin) ->

  F3 = 1.0/3.0
  s = (xin+yin+zin)*F3
  i = Math.floor xin+s
  j = Math.floor yin+s
  k = Math.floor zin+s
  G3 = 1.0/6.0
  t = (i+j+k)*G3
  X0 = i-t
  Y0 = j-t
  Z0 = k-t
  x0 = xin-X0
  y0 = yin-Y0
  z0 = zin-Z0

  if x0 >= y0
    if y0 >= z0
      i1=1; j1=0; k1=0; i2=1; j2=1; k2=0
    else if x0 >= z0
      i1=1; j1=0; k1=0; i2=1; j2=0; k2=1
    else
      i1=0; j1=0; k1=1; i2=1; j2=0; k2=1
  else
    if y0 < z0
      i1=0; j1=0; k1=1; i2=0; j2=1; k2=1
    else if x0 < z0
      i1=0; j1=1; k1=0; i2=0; j2=1; k2=1
    else
      i1=0; j1=1; k1=0; i2=1; j2=1; k2=0

  x1 = x0 - i1 + G3
  y1 = y0 - j1 + G3
  z1 = z0 - k1 + G3
  x2 = x0 - i2 + 2.0*G3
  y2 = y0 - j2 + 2.0*G3
  z2 = z0 - k2 + 2.0*G3
  x3 = x0 - 1.0 + 3.0*G3
  y3 = y0 - 1.0 + 3.0*G3
  z3 = z0 - 1.0 + 3.0*G3

  ii = i & 255
  jj = j & 255
  kk = k & 255

  gi0 = perm[ii+perm[jj+perm[kk]]] % 12
  gi1 = perm[ii+i1+perm[jj+j1+perm[kk+k1]]] % 12
  gi2 = perm[ii+i2+perm[jj+j2+perm[kk+k2]]] % 12
  gi3 = perm[ii+1+perm[jj+1+perm[kk+1]]] % 12

  n0 = n1 = n2 = n3 = 0

  t0 = 0.6 - x0*x0 - y0*y0 - z0*z0
  if t0>0
    t0 *= t0
    n0 = t0 * t0 * dot(x0, y0, z0, grad[gi0])

  t1 = 0.6 - x1*x1 - y1*y1 - z1*z1
  if t1>0
    t1 *= t1
    n1 = t1 * t1 * dot(x1, y1, z1, grad[gi1])

  t2 = 0.6 - x2*x2 - y2*y2 - z2*z2
  if t2>0
    t2 *= t2
    n2 = t2 * t2 * dot(x2, y2, z2, grad[gi2])

  t3 = 0.6 - x3*x3 - y3*y3 - z3*z3
  if t3>0
    t3 *= t3
    n3 = t3 * t3 * dot(x3, y3, z3, grad[gi3])

  16.0*(n0 + n1 + n2 + n3)+1.0

window.simplex_noise = (octaves, x, y, z) ->
  value = 0.0
  for i in [0..octaves-1]
    value += noise(
      x*Math.pow(2, i),
      y*Math.pow(2, i),
      z*Math.pow(2, i))
  value

$ ->
  renderer = new THREE.WebGLRenderer()
  camera = new THREE.PerspectiveCamera 90, 400/300, 0.1, 10000
  scene = new THREE.Scene()
  renderer.setClearColor new THREE.Color(0, 1)
  renderer.setSize 800, 600
  $("body").append(renderer.domElement)

  s = 64
  ts = 512
  $('#canvas').attr({width: ts, height: ts})
  nsh = s/-2

  # Move this in to Sponge? Not sure.
  idx = (x,y,z) ->
    x + y*s + z*s*s

  do ->
    v = new THREE.Vector3(nsh,nsh,nsh)
    THREE.Vector3.prototype.repos = _.partial(THREE.Vector3.prototype.add, v)

  class Sponge
    constructor: (@s) ->
      @length = @s*@s*@s
      @_sponge = new Array(@length)
      @loop (x,y,z) =>
        tv = new THREE.Vector3(x,y,z).multiplyScalar(1/@s)
        tv3 = new THREE.Vector3(tv.x+1,tv.y+1,tv.z+1).multiplyScalar(3)
        tv5 = tv.clone().multiplyScalar(5)

        if tv.y <= 0.8
          plateau_falloff = 1.0
        else if tv.y > 0.8 and tv.y < 0.9
          plateau_falloff = 1.0-(tv.y-0.8)*10.0
        else
          plateau_falloff = 0.0

        center_falloff = 0.1/(Math.pow((tv.x-0.5)*1.5, 2) + Math.pow((tv.y-1.0)*0.8, 2) + Math.pow((tv.z-0.5)*1.5, 2))

        caves = Math.pow simplex_noise(1, tv5.x, tv5.y, tv5.z), 3
        density = simplex_noise(5,tv.x,tv.y/2,tv.z) * plateau_falloff * center_falloff * Math.pow(simplex_noise(1, tv3.x, tv3.y, tv3.z)+0.4, 1.8)
        density = 0 if caves < 0.5
        @_sponge[idx x, y, z] = density > 3.1

    lookup: (x,y,z) ->
      if y == undefined
        @_sponge[x]
      else
        if (x >= 0 and x < @s) and (y >= 0 and y < @s) and (z >= 0 and z < @s) then @_sponge[idx x,y,z] else false

    place: (x,y,z,v) ->
      if z == undefined
        @_sponge[x] = y
      else
      if (x >= 0 and x < @s) and (y >= 0 and y < @s) and (z >= 0 and z < @s) then @_sponge[idx x,y,z] = v else false

    neighbors: (x,y,z) ->
      n = [[x-1,y,z],[x+1,y,z],[x,y-1,z],[x,y+1,z],[x,y,z-1],[x,y,z+1]]
      (@lookup i... for i in n)

    loop: (fn) ->
      for z in [0..@s-1]
        for y in [0..@s-1]
          for x in [0..@s-1]
            fn x, y, z, @lookup(idx(x,y,z))

    loop_flat: (fn) ->
      for i in [0..@length-1]
        fn i, @lookup(i)

  sponge = new Sponge(s)
  console.log sponge

  class Ray
    constructor: (x, y, z) ->
      @vector = new THREE.Vector3 x, y, z
      @rasterize()

    rasterize: ->
      x = @vector.x; y = @vector.y; z = @vector.z
      sx = x*0.2
      sy = y*0.2
      sz = z*0.2
  
      x = if x < 0 then -1 else 1
      y = if y < 0 then -1 else 1
      z = if z < 0 then -1 else 1
  
      cx = 0
      cy = 0
      cz = 0
  
      @offsets = new Array(s)
 
      i = 0 
      while i < s
        ncx = Math.floor x
        ncy = Math.floor y
        ncz = Math.floor z
        if ncx != cx or ncy != cy or ncz != cz
          depth = Math.sqrt x*x+y*y+z*z
          cx = ncx
          cy = ncy
          cz = ncz
          @offsets[i] = new THREE.Vector3(cx, cy, cz)
          i++
        x+=sx
        y+=sy
        z+=sz

  #http://www.softimageblog.com/archives/115
  sphere_dist = (n) ->
    dist = new Array(n)
    inc = Math.PI * (3 - Math.sqrt 5)
    o = 2/n
    for i in [1..n]
      y = i * o - 1 + (o/2)
      r = Math.sqrt(1) - y*y
      phi = i * inc
      dist[i-1] = new Ray(Math.cos(phi) * r, y, Math.sin(phi) * r)
    dist

  rays = sphere_dist 50
#
#  to_canv = (inp)->
#    canvas = document.getElementById('canvas')
#    canvasWidth  = canvas.width
#    canvasHeight = canvas.height
#    ctx = canvas.getContext('2d')
#    imageData = ctx.getImageData(0, 0, canvasWidth, canvasHeight)
#
#    data = imageData.data
#
#    for y in [0..canvasHeight-1]
#      for x in [0..canvasWidth-1]
#        index = (y * canvasWidth + x) * 4
#
#        data[index]   = inp[index]
#        data[++index] = inp[index]
#        data[++index] = inp[index]
#        data[++index] = inp[index]
#
#    ctx.putImageData(imageData, 0, 0)
#
##  readable test data
##  cube = []
##  for x in [0..s-1]
##    for y in [0..s-1]
##      for z in [0..s-1]
##        i = x + y*s + z*s*s
##        cube[i] = [i,i,i]
#
  calc_occs = ->
    r = new THREE.WebGLRenderer({preserveDrawingBuffer: true})
    r.setSize ts, ts
    #$('body').prepend(r.domElement)

    derp = []
    sponge.loop_flat (i,cube) ->
      if cube then derp.push 255,255,255,255 else derp.push 0,0,0,255
    for i in [1..ts*ts*4-derp.length]
      derp.push 0,0,0,255

    input = new Uint8Array(derp)
    console.log "input non zero", _.filter(input,((e,i) -> i%4 != 3 and e > 0)).length
    #to_canv input
    t = new THREE.DataTexture input, ts, ts, THREE.RGBAFormat, THREE.UnsignedByteType, new THREE.UVMapping(), undefined, undefined, THREE.NearestFilter, THREE.NearestFilter
    t.needsUpdate = true

    rs = rays[0].offsets

    sc = new THREE.Scene()
    g = new THREE.PlaneGeometry(2,2)
    m = new THREE.Mesh g, new THREE.ShaderMaterial({
      uniforms: {
        t: {type: 't', value: t},
        s: {type: 'f', value: s},
        ts: {type: 'f', value: ts},
        rayoffs: {type: 'v3v', value: rs}},
      vertexShader: $('#occv').text(),
      fragmentShader: $('#occf').text().replace('rayoffs[64]',"rayoffs[#{s}]").replace('r < 64',"r < #{s}")})
    console.log m
    sc.add m

    #tgt = new THREE.WebGLRenderTarget()
    #tgt.width = tgt.height = s
    output = new Uint8Array ts*ts*4

    occ = new Array(ts*ts)

    to_xyz = (cux,cuy) ->
      fidx = cux + cuy * ts
      coz = Math.floor(fidx/(s*s))
      fidx -= coz * s * s
      coy = Math.floor(fidx/s)
      fidx -= coy * s
      cox = fidx
      return [cox,coy,coz]

    gl = r.getContext()
    for ray,ray_index in rays
      console.log ray_index unless ray_index % 10
      m.material.uniforms.rayoffs.value = ray.offsets

      #r.render sc, camera, tgt
      r.render sc, camera
      gl.readPixels(0,0,ts,ts,gl.RGBA,gl.UNSIGNED_BYTE,output)

      # validate coordinate xforms
      #sponge.loop (x,y,z,i) ->
      #  [x1,y1,z1] = to_xyz i % ts, Math.floor(i/ts)
      #  console.log "fail x", i unless x == x1 and y == y1 and z == z1

      # validate texture read
      #failed = 0
      #for op,opidx in output
      #  cidx = Math.floor opidx/4
      #  xyz = to_xyz(cidx % ts, Math.floor(cidx/ts))
      #  unless op == input[opidx]
      #    #console.log "fail", opidx, op, input[opidx], xyz, (i for i in xyz)
      #    console.log "fail", opidx, "out:", op, "in:", input[opidx]
      #    failed++
      #  #else
      #  #  console.log "good", opidx
      #console.log "total failed", failed
      #to_canv output
      for comp,l in output by 4
        l4 = l/4
        occ[l4] = [[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]] if ray_index == 0
        x = ray.vector.x; y = ray.vector.y; z = ray.vector.z
        if x < 0
          occ[l4][0][0] -= x
          occ[l4][0][1] -= x if comp > 0
        if x > 0
          occ[l4][1][0] += x
          occ[l4][1][1] += x if comp > 0
        if y < 0
          occ[l4][2][0] -= y
          occ[l4][2][1] -= y if comp > 0
        if y > 0
          occ[l4][3][0] += y
          occ[l4][3][1] += y if comp > 0
        if z < 0
          occ[l4][4][0] -= z
          occ[l4][4][1] -= z if comp > 0
        if z > 0
          occ[l4][5][0] += z
          occ[l4][5][1] += z if comp > 0

    o = []
    for cube_faces in occ
      t = []
      for face in cube_faces
        t.push if face[0] and face[1] then face[1]/face[0] else 0
      o.push t
    o

  occs_per_node = calc_occs()

  g = new THREE.Geometry()

  idxvs = (x,y,z) ->
    g.vertices.push new THREE.Vector3(x,y,z).repos()
    g.vertices.length - 1

  occs_per_vert = []
  sponge.loop (x,y,z,cube) ->
    if cube
      _.each sponge.neighbors(x,y,z), (e2,i2) -> 
        unless e2[0]
          f = switch i2
            when 0 then new THREE.Face4 idxvs(x,y,z), idxvs(x,y,z+1), idxvs(x,y+1,z+1), idxvs(x,y+1,z)
            when 1 then new THREE.Face4 idxvs(x+1,y,z), idxvs(x+1,y+1,z), idxvs(x+1,y+1,z+1), idxvs(x+1,y,z+1)
            when 2 then new THREE.Face4 idxvs(x,y,z), idxvs(x+1,y,z), idxvs(x+1,y,z+1), idxvs(x,y,z+1)
            when 3 then new THREE.Face4 idxvs(x,y+1,z), idxvs(x,y+1,z+1), idxvs(x+1,y+1,z+1), idxvs(x+1,y+1,z)
            when 4 then new THREE.Face4 idxvs(x,y,z), idxvs(x,y+1,z), idxvs(x+1,y+1,z), idxvs(x+1,y,z)
            when 5 then new THREE.Face4 idxvs(x,y,z+1), idxvs(x+1,y,z+1), idxvs(x+1,y+1,z+1), idxvs(x,y+1,z+1)
          g.faces.push f
          g.faceVertexUvs[0].push [new THREE.Vector2(0,0),new THREE.Vector2(0,1),new THREE.Vector2(1,1),new THREE.Vector2(1,0)]
          for _ in [1..4]
            occs_per_vert.push occs_per_node[idx x,y,z][i2]

  console.log g if s < 5

  m = new THREE.Mesh g, new THREE.ShaderMaterial({
    vertexShader: $('#shlightv').text(),
    fragmentShader: $('#shlightf').text(),
    uniforms: {t: {type: 't', value: THREE.ImageUtils.loadTexture 'stone.jpg'}, ScaleFactor: 1.0},
    attributes: {ocv: {type: 'f', value: occs_per_vert}}})
  scene.add m

  g.computeFaceNormals()
  g.computeVertexNormals()
  g.normalsNeedUpdate = true
  g.uvsNeedUpdate = true
  console.log g

  xx = yy = zz = 0
  $('body').keydown (e) ->
    switch e.which
      when 87 then zz = -1
      when 83 then zz = 1
      when 65 then xx = -1
      when 68 then xx = 1
      when 32 then yy = 1
  $('body').keyup (e) ->
    switch e.which
      when 87 then zz = 0
      when 83 then zz = 0
      when 65 then xx = 0
      when 68 then xx = 0
      when 32 then yy = 0
  track = false
  last = pos = [0,0]
  dx = dy = 0
  $('body').mousedown (e) ->
    track = true;
  $('body').mouseup (e) ->
    track = false;
  $('body').mousemove (e) ->
    last = pos
    pos = [e.screenX,e.screenY]

  samples = 0
  update_deltas = ->
    samples++
    if track
      return unless samples % 10
      dx = (last[0] - pos[0]) / (last[0] - pos[0] * 3)
      dy = (last[1] - pos[1]) / (last[1] - pos[1] * 3)
    else
      dx = dy = 0

  camera.position.z = s * 0.75
  camera.eulerOrder = 'YXZ'
  precalc = ([s * Math.sin(i), s * Math.cos(i)] for i in [0.0..2*Math.PI] by 0.0025)
  doit = (t) ->
    u = precalc[Math.floor(t % precalc.length)] or [0,1]
    camera.translateX xx if xx
    camera.translateY yy if yy
    camera.translateZ zz if zz
    update_deltas()
    camera.rotation.y += dx
    camera.rotation.x += dy
    renderer.render scene, camera
    requestAnimationFrame doit

  doit 0
