Loaded 260 links
# UCSC hub corrections (replaces earlier file listing)

Three real errors in the previous "UCSC Exact Files" comment, all caught by review:

1. **`type chain` is not a valid track-hub type.** Track hubs need indexed binary `bigChain` files, not gzipped text. Conversion done — every chain has been converted to a `.bigChain.bb` + `.bigChain.link.bb` pair (UCSC bigBed indexed) via the standard `bigChain.as` / `bigLink.as` schema. The raw `.chain.gz` files remain on Dropbox + Git as download-only artifacts.

2. **`genomes.txt` was missing real `defaultPos` and the `.2bit` files were not in place.** Every hub now jumps to its **dhps locus** (PVP01_1429500 ortholog on chr14) on first load, and every hub directory has a `.2bit` for its assembly.

3. **bigMaf cannot share a composite with chain tracks** — UCSC composites must be a single `type`. The multi-z is now a standalone track at the top of each `trackDb.txt`; chains live under a separate `brc_pangenome_chains` composite of type `bigChain`.

Per-hub `trackDb.txt` now has the shape:

```
track {name}_multiz                          # standalone, type bigMaf
track brc_pangenome_chains                   # composite, type bigChain  (7 sub-tracks)
track brc_pangenome_annot                    # composite, type bigBed 12 (4 sub-tracks)
track brc_pangenome_select                   # composite, PvP01 only     (3 sub-tracks)
```

## defaultPos per hub

| Hub | dhps region (chr14 ortholog) |
|---|---|
| PvP01 (`GCA_900093555.2`) | `LT635625.2:1264700-1277700` |
| Sal-I (`GCA_000002415.2`) | `CM000455.1:1251700-1264600` |
| PvW1 (`GCA_914969965.1`) | `CAJZCX010000003.1:1325900-1338900` |
| PAM (`GCA_949152365.1`) | `CASCJQ010000014.1:1380400-1393500` |
| PvSY56 (`GCA_003402215.1`) | `QMFC01000014.1:1260000-1273000` |
| PvT01 (`GCA_900093545.1`) | `LT615252.1:1275900-1288800` |
| PvC01 (`GCA_900093535.1`) | `LT615269.1:1287800-1300700` |
| MHC087 (`GCA_040114635.1`) | `JBDKXN010000002.1:1285000-1298000` |

## New / changed files on Dropbox

- **Top-level**: [genomes.txt](https://www.dropbox.com/scl/fi/g4g1k1ky2utdraa7i4fzp/genomes.txt?rlkey=ja58vv6hki7lqsh19xr0rnk9q&dl=0) ← rewritten

**Per hub** (8 of them):

| File | Was | Now |
|---|---|---|
| `trackDb.txt` | mixed-type composite, `type chain` | standalone bigMaf + `bigChain` composite |
| `{ACC}.2bit` | missing | symlinked from `projection/A2_kegalign/2bit/` |
| `chains/{ACC}_to_{TGT}.bigChain.bb` | did not exist | bigBed 6+6 (10–400 KB each) |
| `chains/{ACC}_to_{TGT}.bigChain.link.bb` | did not exist | bigBed 4+1 link companion (600 KB – 1.8 MB each) |
| `chains/{ACC}_to_{TGT}.chain.gz` | a track-hub track (broken) | raw chain, download-only |

## Per-hub file listing (clickable)

### PvP01 — `GCA_900093555.2/`

| File | Link |
|---|---|
| `trackDb.txt` | [Dropbox](https://www.dropbox.com/scl/fi/nvvjf7iw6v87e2uftjvzf/trackDb.txt?rlkey=t6vbnvui4ouwuqvt13rt4woe2&dl=0) |
| `GCA_900093555.2.2bit` | [Dropbox](https://www.dropbox.com/scl/fi/jotwcgoyfw6klwa4o7r4l/GCA_900093555.2.2bit?rlkey=qg7huerh8db5n6s2egd387bre&dl=0) |
| `GCA_900093555.2_to_GCA_000002415.2.bigChain.bb` (→ Sal-I) | [bb](https://www.dropbox.com/scl/fi/cpas475bq2ln7mw0z3ae7/GCA_900093555.2_to_GCA_000002415.2.bigChain.bb?rlkey=w85xzjexjffghqf0dsg6xmuli&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/5q45p2h3jxi170nmlvl0w/GCA_900093555.2_to_GCA_000002415.2.bigChain.link.bb?rlkey=hzhuyi87rq1noq21w3t7cd325&dl=0) |
| `GCA_900093555.2_to_GCA_914969965.1.bigChain.bb` (→ PvW1) | [bb](https://www.dropbox.com/scl/fi/pyb62hl8awsw8sk6sy301/GCA_900093555.2_to_GCA_914969965.1.bigChain.bb?rlkey=30dxidf2fk8d47jritrejuek7&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/ps7eqp0iskcfppku8c648/GCA_900093555.2_to_GCA_914969965.1.bigChain.link.bb?rlkey=njy87e7zil4xpcnf8ftjf9pgg&dl=0) |
| `GCA_900093555.2_to_GCA_949152365.1.bigChain.bb` (→ PAM) | [bb](https://www.dropbox.com/scl/fi/sky57z5za7yy57asabjtt/GCA_900093555.2_to_GCA_949152365.1.bigChain.bb?rlkey=ta7qn0v1pq218nt1s393v2jdv&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/a3bjf27jsk589qiuknwqp/GCA_900093555.2_to_GCA_949152365.1.bigChain.link.bb?rlkey=3w0qqiioq32slnbva2yj4hykf&dl=0) |
| `GCA_900093555.2_to_GCA_003402215.1.bigChain.bb` (→ PvSY56) | [bb](https://www.dropbox.com/scl/fi/nl0jb1v8inkalgaukq1m9/GCA_900093555.2_to_GCA_003402215.1.bigChain.bb?rlkey=456jtfgvm18qeyb9nirj50xmk&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/b8j4t58pg4m82yktwl8a0/GCA_900093555.2_to_GCA_003402215.1.bigChain.link.bb?rlkey=eqt3l1ds9s8ac37sstyqkzctg&dl=0) |
| `GCA_900093555.2_to_GCA_900093545.1.bigChain.bb` (→ PvT01) | [bb](https://www.dropbox.com/scl/fi/94h2oglvuop9t2dhb18gl/GCA_900093555.2_to_GCA_900093545.1.bigChain.bb?rlkey=3ew96jssdw25xr7575d26kkv1&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/8rfpmhsvj5pilxuwfh35u/GCA_900093555.2_to_GCA_900093545.1.bigChain.link.bb?rlkey=xsp6rt5z4zjrdny4lf30ltrhj&dl=0) |
| `GCA_900093555.2_to_GCA_900093535.1.bigChain.bb` (→ PvC01) | [bb](https://www.dropbox.com/scl/fi/mrvqtcj9vgy45wlz3oocx/GCA_900093555.2_to_GCA_900093535.1.bigChain.bb?rlkey=sn68wpf17qv4u83v019sepl5z&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/ucgh34b5fdkgnhjxkbogu/GCA_900093555.2_to_GCA_900093535.1.bigChain.link.bb?rlkey=arwe1dn93swx4ut92aer5bzys&dl=0) |
| `GCA_900093555.2_to_GCA_040114635.1.bigChain.bb` (→ MHC087) | [bb](https://www.dropbox.com/scl/fi/e5grai7njoxv4jbbum86b/GCA_900093555.2_to_GCA_040114635.1.bigChain.bb?rlkey=f4fvaluezv37wdpg1nr4nfsb4&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/ws7e448sjxfd00vux7x3g/GCA_900093555.2_to_GCA_040114635.1.bigChain.link.bb?rlkey=7em712fea9o0fpq2db75yq99z&dl=0) |

### Sal-I — `GCA_000002415.2/`

| File | Link |
|---|---|
| `trackDb.txt` | [Dropbox](https://www.dropbox.com/scl/fi/4lr4tu19hmncsspag29kc/trackDb.txt?rlkey=ittvke5grd8166k687vej8zae&dl=0) |
| `GCA_000002415.2.2bit` | [Dropbox](https://www.dropbox.com/scl/fi/0k4x2zesnfeqvyij0x274/GCA_000002415.2.2bit?rlkey=berwe02jbba21641q8sn0fx7i&dl=0) |
| `GCA_000002415.2_to_GCA_900093555.2.bigChain.bb` (→ PvP01) | [bb](https://www.dropbox.com/scl/fi/uy99yjuajopb4o2yxf3ga/GCA_000002415.2_to_GCA_900093555.2.bigChain.bb?rlkey=dqhmnx3qp2mzwxuwn6jrtk525&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/79mqtxjfvsu1w1zrlqt7o/GCA_000002415.2_to_GCA_900093555.2.bigChain.link.bb?rlkey=vtb0jxxvbxr9fj3x9je5pgkgq&dl=0) |
| `GCA_000002415.2_to_GCA_914969965.1.bigChain.bb` (→ PvW1) | [bb](https://www.dropbox.com/scl/fi/qdj6k5oz89rh72shxulc4/GCA_000002415.2_to_GCA_914969965.1.bigChain.bb?rlkey=q95v7e8qpru0lj41x1uq08aub&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/dug2v27foftnsh6l1gy2y/GCA_000002415.2_to_GCA_914969965.1.bigChain.link.bb?rlkey=6j5n681310xuac93zg1rf7n7m&dl=0) |
| `GCA_000002415.2_to_GCA_949152365.1.bigChain.bb` (→ PAM) | [bb](https://www.dropbox.com/scl/fi/ik2ihk7wrmr2blt424io5/GCA_000002415.2_to_GCA_949152365.1.bigChain.bb?rlkey=mzvwnteowst4ksq8grst28prq&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/3uovi6b6p5so2zps6oj1l/GCA_000002415.2_to_GCA_949152365.1.bigChain.link.bb?rlkey=6ngot9himpevj48y160ii0ws0&dl=0) |
| `GCA_000002415.2_to_GCA_003402215.1.bigChain.bb` (→ PvSY56) | [bb](https://www.dropbox.com/scl/fi/3bmkzlut7jakq828501n7/GCA_000002415.2_to_GCA_003402215.1.bigChain.bb?rlkey=dh5icr5olu3aemycdgep3e2ud&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/fnevggpxgsakcn9kgrlz0/GCA_000002415.2_to_GCA_003402215.1.bigChain.link.bb?rlkey=d3uruy76hkwk8kzqm51fi9b7c&dl=0) |
| `GCA_000002415.2_to_GCA_900093545.1.bigChain.bb` (→ PvT01) | [bb](https://www.dropbox.com/scl/fi/u1ve85qr9lm78if75ov1n/GCA_000002415.2_to_GCA_900093545.1.bigChain.bb?rlkey=zaxqv9nqgnh7gd2tkfr7ymegk&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/f35wymelzygcfh4ttsm6r/GCA_000002415.2_to_GCA_900093545.1.bigChain.link.bb?rlkey=09i0sgu1fokiocg89z673n9u1&dl=0) |
| `GCA_000002415.2_to_GCA_900093535.1.bigChain.bb` (→ PvC01) | [bb](https://www.dropbox.com/scl/fi/4buzycoy54z1n9q95zsvi/GCA_000002415.2_to_GCA_900093535.1.bigChain.bb?rlkey=adxonjmnhpsy0482pq4oyjkg8&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/a7nrmg6zw9wbvfir5kcf4/GCA_000002415.2_to_GCA_900093535.1.bigChain.link.bb?rlkey=buziplelv4l3ew4k3mnaru5wc&dl=0) |
| `GCA_000002415.2_to_GCA_040114635.1.bigChain.bb` (→ MHC087) | [bb](https://www.dropbox.com/scl/fi/mfzhbc47rnqqsw4utaw2x/GCA_000002415.2_to_GCA_040114635.1.bigChain.bb?rlkey=wp81b04u0whyi0impr8qlt13o&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/lyp14t32dfdvovlwuofw8/GCA_000002415.2_to_GCA_040114635.1.bigChain.link.bb?rlkey=hq67x0vhkj1jl63cenh6nqbc1&dl=0) |

### PvW1 — `GCA_914969965.1/`

| File | Link |
|---|---|
| `trackDb.txt` | [Dropbox](https://www.dropbox.com/scl/fi/k7zgr53ejqdeaifsv86aw/trackDb.txt?rlkey=irckckz81o2e4n3f2kqjg1nio&dl=0) |
| `GCA_914969965.1.2bit` | [Dropbox](https://www.dropbox.com/scl/fi/ljzux8xufta0km82e19s8/GCA_914969965.1.2bit?rlkey=ds8i0hkmtkr5xa4228b6pk75d&dl=0) |
| `GCA_914969965.1_to_GCA_900093555.2.bigChain.bb` (→ PvP01) | [bb](https://www.dropbox.com/scl/fi/mizcr987fzd5goghx9vwd/GCA_914969965.1_to_GCA_900093555.2.bigChain.bb?rlkey=5s5dxyyhoq6q97hme071p1hbf&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/m2d342ze5oisx2tfgshe8/GCA_914969965.1_to_GCA_900093555.2.bigChain.link.bb?rlkey=n66sn2205iq5t7ubv8cqrnqel&dl=0) |
| `GCA_914969965.1_to_GCA_000002415.2.bigChain.bb` (→ Sal-I) | [bb](https://www.dropbox.com/scl/fi/ce3n72xijetptyjdqwel8/GCA_914969965.1_to_GCA_000002415.2.bigChain.bb?rlkey=0sd556f0uwuqvibbek0jp9w3i&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/v8t9uilpedwphxin8qo0p/GCA_914969965.1_to_GCA_000002415.2.bigChain.link.bb?rlkey=jvvtfpvcr2h0c4q0ns6ohfh3p&dl=0) |
| `GCA_914969965.1_to_GCA_949152365.1.bigChain.bb` (→ PAM) | [bb](https://www.dropbox.com/scl/fi/88oqr8d9f45aoayflebq3/GCA_914969965.1_to_GCA_949152365.1.bigChain.bb?rlkey=33d5jyejroot241zys0eyxn4b&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/us9eetpz5pdf0jy71t1pg/GCA_914969965.1_to_GCA_949152365.1.bigChain.link.bb?rlkey=wtm9qt0fwlimlecvspcwkyj9a&dl=0) |
| `GCA_914969965.1_to_GCA_003402215.1.bigChain.bb` (→ PvSY56) | [bb](https://www.dropbox.com/scl/fi/723gbyh27po8zxttpytdv/GCA_914969965.1_to_GCA_003402215.1.bigChain.bb?rlkey=nzih80t1zy6sdmxb4m8uiqgl6&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/nzm7kk7x16qdsxk6r9ol2/GCA_914969965.1_to_GCA_003402215.1.bigChain.link.bb?rlkey=x0j63vwkifpybwxlsad75ic13&dl=0) |
| `GCA_914969965.1_to_GCA_900093545.1.bigChain.bb` (→ PvT01) | [bb](https://www.dropbox.com/scl/fi/oqsvkar76gt2qmazi2ug9/GCA_914969965.1_to_GCA_900093545.1.bigChain.bb?rlkey=xq6sxg265vcdhrxlfl25vu871&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/i8d3de7dfi8ywzi93wxki/GCA_914969965.1_to_GCA_900093545.1.bigChain.link.bb?rlkey=cb95n45wfqk0k32je6tn61t5g&dl=0) |
| `GCA_914969965.1_to_GCA_900093535.1.bigChain.bb` (→ PvC01) | [bb](https://www.dropbox.com/scl/fi/cqbnvdguh59fy3va3vzrb/GCA_914969965.1_to_GCA_900093535.1.bigChain.bb?rlkey=k0kxvr0ofdpz61ducv00gxczv&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/ra9ov8d2k3zs71zqu5l8f/GCA_914969965.1_to_GCA_900093535.1.bigChain.link.bb?rlkey=k5v7k42uvzsod7rgwp9b2cgrg&dl=0) |
| `GCA_914969965.1_to_GCA_040114635.1.bigChain.bb` (→ MHC087) | [bb](https://www.dropbox.com/scl/fi/7ay7v2c3varuc5yxl7dpy/GCA_914969965.1_to_GCA_040114635.1.bigChain.bb?rlkey=ktxzd6da9srtiwj19jnujj8zn&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/ofiwpw78d2v17240wh0qk/GCA_914969965.1_to_GCA_040114635.1.bigChain.link.bb?rlkey=yejb9lshm7suvtfhiiea1yiqv&dl=0) |

### PAM — `GCA_949152365.1/`

| File | Link |
|---|---|
| `trackDb.txt` | [Dropbox](https://www.dropbox.com/scl/fi/glr8icfwlf38zcrh1jcge/trackDb.txt?rlkey=lrh1bys9bh2nag90xofp934kd&dl=0) |
| `GCA_949152365.1.2bit` | [Dropbox](https://www.dropbox.com/scl/fi/138c1iunqm5cmdf3oqgv9/GCA_949152365.1.2bit?rlkey=92du1e9l0fual3n7gm1tmfvbw&dl=0) |
| `GCA_949152365.1_to_GCA_900093555.2.bigChain.bb` (→ PvP01) | [bb](https://www.dropbox.com/scl/fi/asplozl97c4mlosczxqs3/GCA_949152365.1_to_GCA_900093555.2.bigChain.bb?rlkey=uzmjgavzo98acj6a1whqfxbed&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/gnakkykpushheennxujz5/GCA_949152365.1_to_GCA_900093555.2.bigChain.link.bb?rlkey=1jy6o1wis3it24hqucfkpw45w&dl=0) |
| `GCA_949152365.1_to_GCA_000002415.2.bigChain.bb` (→ Sal-I) | [bb](https://www.dropbox.com/scl/fi/6bd2hik713zt3i8wcqtxf/GCA_949152365.1_to_GCA_000002415.2.bigChain.bb?rlkey=ibt31p79rlciq7qxx9i9jc99m&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/snjg2ddibzh7salx976c6/GCA_949152365.1_to_GCA_000002415.2.bigChain.link.bb?rlkey=r1nzyby5ss0at4qc62p67xqn6&dl=0) |
| `GCA_949152365.1_to_GCA_914969965.1.bigChain.bb` (→ PvW1) | [bb](https://www.dropbox.com/scl/fi/aao6h0qnfuep9ccgtwc95/GCA_949152365.1_to_GCA_914969965.1.bigChain.bb?rlkey=n1p9spelhdbs6vihdtq3qwtaw&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/jkaonms38d5zyicx8lsti/GCA_949152365.1_to_GCA_914969965.1.bigChain.link.bb?rlkey=2ij0siudfr8qhavrajww4eoem&dl=0) |
| `GCA_949152365.1_to_GCA_003402215.1.bigChain.bb` (→ PvSY56) | [bb](https://www.dropbox.com/scl/fi/xftcalrqar3pziqtidcgc/GCA_949152365.1_to_GCA_003402215.1.bigChain.bb?rlkey=uf690lb2pssfa04kvjkh25jou&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/l316noep9hyg6jjotmqby/GCA_949152365.1_to_GCA_003402215.1.bigChain.link.bb?rlkey=jhcya5jpajeo1t7a4r07hkc66&dl=0) |
| `GCA_949152365.1_to_GCA_900093545.1.bigChain.bb` (→ PvT01) | [bb](https://www.dropbox.com/scl/fi/hyjch60qy5lb79ojfoan3/GCA_949152365.1_to_GCA_900093545.1.bigChain.bb?rlkey=2etq9mir2qye87u5ey0ojv201&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/g0efs3t18rc610qm8qq9o/GCA_949152365.1_to_GCA_900093545.1.bigChain.link.bb?rlkey=5qy27m0145yqiooif8tn40jff&dl=0) |
| `GCA_949152365.1_to_GCA_900093535.1.bigChain.bb` (→ PvC01) | [bb](https://www.dropbox.com/scl/fi/8rlzn1ob0d4ldobp6165l/GCA_949152365.1_to_GCA_900093535.1.bigChain.bb?rlkey=k49pyi578kxsdmhidxn2r7lsl&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/qgyuu0dhl58oh24yg802r/GCA_949152365.1_to_GCA_900093535.1.bigChain.link.bb?rlkey=wwewa3gc3cmqm9esdc912z3h1&dl=0) |
| `GCA_949152365.1_to_GCA_040114635.1.bigChain.bb` (→ MHC087) | [bb](https://www.dropbox.com/scl/fi/otf3535umna32755mmnqs/GCA_949152365.1_to_GCA_040114635.1.bigChain.bb?rlkey=qpkyd5dwm014ae2nime9sbz78&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/uu6uzect0qwpg1o26ncl8/GCA_949152365.1_to_GCA_040114635.1.bigChain.link.bb?rlkey=zwxscqvj6etwgg3w2375aftyu&dl=0) |

### PvSY56 — `GCA_003402215.1/`

| File | Link |
|---|---|
| `trackDb.txt` | [Dropbox](https://www.dropbox.com/scl/fi/ljnvcuaicenbju4uhu3o0/trackDb.txt?rlkey=6zem840o55lryn1mqe4nj3439&dl=0) |
| `GCA_003402215.1.2bit` | [Dropbox](https://www.dropbox.com/scl/fi/tep53ardq6zdtioys0g7h/GCA_003402215.1.2bit?rlkey=gcchfvim45t32s3fcbywe9acp&dl=0) |
| `GCA_003402215.1_to_GCA_900093555.2.bigChain.bb` (→ PvP01) | [bb](https://www.dropbox.com/scl/fi/oz72k3r7ftq88xbqc28ct/GCA_003402215.1_to_GCA_900093555.2.bigChain.bb?rlkey=ot8arpn9f4jez4srhwybznpkr&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/hmmexl78pf7i3altkpa61/GCA_003402215.1_to_GCA_900093555.2.bigChain.link.bb?rlkey=xncvfnmp1qx2989k8ucghi80d&dl=0) |
| `GCA_003402215.1_to_GCA_000002415.2.bigChain.bb` (→ Sal-I) | [bb](https://www.dropbox.com/scl/fi/6iman9voeamr815prf6ui/GCA_003402215.1_to_GCA_000002415.2.bigChain.bb?rlkey=4ctd92pafmg3m7ytw9ijlfbgv&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/dkwj9zamwn86g7a4ygyz2/GCA_003402215.1_to_GCA_000002415.2.bigChain.link.bb?rlkey=vaxmyndi31rnpp2cpjhjjlw94&dl=0) |
| `GCA_003402215.1_to_GCA_914969965.1.bigChain.bb` (→ PvW1) | [bb](https://www.dropbox.com/scl/fi/2s9hcxwa774lwh4c04seh/GCA_003402215.1_to_GCA_914969965.1.bigChain.bb?rlkey=n65a3fuaq3solb9g2c5mdteq2&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/mmx4dvwqnzovpbdibksea/GCA_003402215.1_to_GCA_914969965.1.bigChain.link.bb?rlkey=qol6j6vp85qlv0l0hntbxbts9&dl=0) |
| `GCA_003402215.1_to_GCA_949152365.1.bigChain.bb` (→ PAM) | [bb](https://www.dropbox.com/scl/fi/2vsfx2dfqrdyg092bcyht/GCA_003402215.1_to_GCA_949152365.1.bigChain.bb?rlkey=5u5ke1oof19eyanxn2f1jobmt&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/bmatkadiawkeqjbg7iogi/GCA_003402215.1_to_GCA_949152365.1.bigChain.link.bb?rlkey=1jqikik5qxticerpenhkqnizk&dl=0) |
| `GCA_003402215.1_to_GCA_900093545.1.bigChain.bb` (→ PvT01) | [bb](https://www.dropbox.com/scl/fi/edkbatcn4ueeddcg627qq/GCA_003402215.1_to_GCA_900093545.1.bigChain.bb?rlkey=kb7vwizmxk5wtpy0mc66k49v4&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/0nsva7llxv31frtouxa89/GCA_003402215.1_to_GCA_900093545.1.bigChain.link.bb?rlkey=vukk9637yr0k5rxl9nku8zjk2&dl=0) |
| `GCA_003402215.1_to_GCA_900093535.1.bigChain.bb` (→ PvC01) | [bb](https://www.dropbox.com/scl/fi/00vzynm1kccuzpp978tym/GCA_003402215.1_to_GCA_900093535.1.bigChain.bb?rlkey=yy59mb99esqt541qkty98e9of&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/sg786sa1twfctxiracu4y/GCA_003402215.1_to_GCA_900093535.1.bigChain.link.bb?rlkey=17chse16sm1nxuhzfbj3ivx7v&dl=0) |
| `GCA_003402215.1_to_GCA_040114635.1.bigChain.bb` (→ MHC087) | [bb](https://www.dropbox.com/scl/fi/stxhntk0gcotqh0l7imug/GCA_003402215.1_to_GCA_040114635.1.bigChain.bb?rlkey=k4727755yycuxp9atisgejs03&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/ymo6jxl9fukhdcsfr5asz/GCA_003402215.1_to_GCA_040114635.1.bigChain.link.bb?rlkey=rbls4bhd1g2woj4i219b315y2&dl=0) |

### PvT01 — `GCA_900093545.1/`

| File | Link |
|---|---|
| `trackDb.txt` | [Dropbox](https://www.dropbox.com/scl/fi/wgrjd7e2a3vur5eege1ot/trackDb.txt?rlkey=11n8yevq4f1b4t5ryo0ioo6vq&dl=0) |
| `GCA_900093545.1.2bit` | [Dropbox](https://www.dropbox.com/scl/fi/lbgb0mimcw1zfj9cidicc/GCA_900093545.1.2bit?rlkey=h6tk19zrck3zh3e8cj6533nv6&dl=0) |
| `GCA_900093545.1_to_GCA_900093555.2.bigChain.bb` (→ PvP01) | [bb](https://www.dropbox.com/scl/fi/xextdinqkb1m7k6kbw8vk/GCA_900093545.1_to_GCA_900093555.2.bigChain.bb?rlkey=uewtgvbztqjxjopyft0vp7x4z&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/zb1cddgktvxqrhzwgsifs/GCA_900093545.1_to_GCA_900093555.2.bigChain.link.bb?rlkey=y5bsgz29hl5g4h1ld38g8vn7s&dl=0) |
| `GCA_900093545.1_to_GCA_000002415.2.bigChain.bb` (→ Sal-I) | [bb](https://www.dropbox.com/scl/fi/9xw5cm4zeweumd7xz9nfb/GCA_900093545.1_to_GCA_000002415.2.bigChain.bb?rlkey=rywjaj99zvtrxk0wbpagru8y8&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/p87c1vdjywlgmuvj94wtq/GCA_900093545.1_to_GCA_000002415.2.bigChain.link.bb?rlkey=t3vdg5kicfzjvbcka0zadeexk&dl=0) |
| `GCA_900093545.1_to_GCA_914969965.1.bigChain.bb` (→ PvW1) | [bb](https://www.dropbox.com/scl/fi/p3wa5dek2f71o08rphcuu/GCA_900093545.1_to_GCA_914969965.1.bigChain.bb?rlkey=gg84fj4dp2c4jmakvq6n37jzj&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/s139xbm91ch6kuqyt77gn/GCA_900093545.1_to_GCA_914969965.1.bigChain.link.bb?rlkey=ogjgl9ngatazzh4nv9hy5erif&dl=0) |
| `GCA_900093545.1_to_GCA_949152365.1.bigChain.bb` (→ PAM) | [bb](https://www.dropbox.com/scl/fi/bakc44x5cwl4j5y6ovy9q/GCA_900093545.1_to_GCA_949152365.1.bigChain.bb?rlkey=jv4ggceiqy81s6u2po8gypjf0&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/z9xng5pubnnbsh1trhusm/GCA_900093545.1_to_GCA_949152365.1.bigChain.link.bb?rlkey=c1ha3k8085tppb75hjagrfefh&dl=0) |
| `GCA_900093545.1_to_GCA_003402215.1.bigChain.bb` (→ PvSY56) | [bb](https://www.dropbox.com/scl/fi/jyf7068fx0jyv4m9tbbfn/GCA_900093545.1_to_GCA_003402215.1.bigChain.bb?rlkey=06vw88ojiyv4lww5qs38kvvuk&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/y5qf7p9gbhcndmg3ctcbi/GCA_900093545.1_to_GCA_003402215.1.bigChain.link.bb?rlkey=pu8xevvxvn15hc4phh6lx73er&dl=0) |
| `GCA_900093545.1_to_GCA_900093535.1.bigChain.bb` (→ PvC01) | [bb](https://www.dropbox.com/scl/fi/rg5ov6looibfp6ew196f8/GCA_900093545.1_to_GCA_900093535.1.bigChain.bb?rlkey=5w49wyi8jdophjyaztvfdx9eb&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/m0qtehozxgh3b75ajnis8/GCA_900093545.1_to_GCA_900093535.1.bigChain.link.bb?rlkey=7w6ya7fmfmcw9dmzozktza754&dl=0) |
| `GCA_900093545.1_to_GCA_040114635.1.bigChain.bb` (→ MHC087) | [bb](https://www.dropbox.com/scl/fi/fqrnu9woclbf5pes2bulq/GCA_900093545.1_to_GCA_040114635.1.bigChain.bb?rlkey=7anycqgcp6iulmry0qlzecsxy&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/uti7kynbqxfdc2sl27f8j/GCA_900093545.1_to_GCA_040114635.1.bigChain.link.bb?rlkey=j03afwuo374y733ujnku7208p&dl=0) |

### PvC01 — `GCA_900093535.1/`

| File | Link |
|---|---|
| `trackDb.txt` | [Dropbox](https://www.dropbox.com/scl/fi/cml6coefnoqcz2vnln1wi/trackDb.txt?rlkey=6b3k6pv9j0duxutw1yeorjuxe&dl=0) |
| `GCA_900093535.1.2bit` | [Dropbox](https://www.dropbox.com/scl/fi/0823kystmbh8ko4x0g3h1/GCA_900093535.1.2bit?rlkey=4n2mcgu1tt9yvnzxuqzrkmpyk&dl=0) |
| `GCA_900093535.1_to_GCA_900093555.2.bigChain.bb` (→ PvP01) | [bb](https://www.dropbox.com/scl/fi/p0f8aj4823mgueq805uri/GCA_900093535.1_to_GCA_900093555.2.bigChain.bb?rlkey=2bgs9bsis385necip7gz3titi&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/0ywswjd9l6epjrwzjfqjw/GCA_900093535.1_to_GCA_900093555.2.bigChain.link.bb?rlkey=jhtix7cdqw2lk937107yvkiut&dl=0) |
| `GCA_900093535.1_to_GCA_000002415.2.bigChain.bb` (→ Sal-I) | [bb](https://www.dropbox.com/scl/fi/ynxpsn7x2ainb61ol1zdt/GCA_900093535.1_to_GCA_000002415.2.bigChain.bb?rlkey=knmtsfvumdxxtsj09i2arp6bi&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/i6ls4isml4usiqdwrvu6z/GCA_900093535.1_to_GCA_000002415.2.bigChain.link.bb?rlkey=qgyct4i3ndd4pwdyk1kpe434k&dl=0) |
| `GCA_900093535.1_to_GCA_914969965.1.bigChain.bb` (→ PvW1) | [bb](https://www.dropbox.com/scl/fi/z0z6xma7cwdvwy4i4e2s9/GCA_900093535.1_to_GCA_914969965.1.bigChain.bb?rlkey=qzc5s4uvtd1gbpcae4m5wtmfa&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/rt1mbaj6jmrk6wgj0hca6/GCA_900093535.1_to_GCA_914969965.1.bigChain.link.bb?rlkey=4mtspddbu9rvqwbln73n8otbs&dl=0) |
| `GCA_900093535.1_to_GCA_949152365.1.bigChain.bb` (→ PAM) | [bb](https://www.dropbox.com/scl/fi/5t6s6mmq1fx9lrfgdqbqp/GCA_900093535.1_to_GCA_949152365.1.bigChain.bb?rlkey=2frjbkx0aq6alksjnie2exmjt&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/zuhexrsmx88nh1ffobas4/GCA_900093535.1_to_GCA_949152365.1.bigChain.link.bb?rlkey=cw0v4yynzdtgp8g934zf68fck&dl=0) |
| `GCA_900093535.1_to_GCA_003402215.1.bigChain.bb` (→ PvSY56) | [bb](https://www.dropbox.com/scl/fi/qw5vmnv84qz78kgylmahq/GCA_900093535.1_to_GCA_003402215.1.bigChain.bb?rlkey=e1vevktmbrv7buut5d3s5c10m&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/ho6ibhxtjbiro9xapcb8o/GCA_900093535.1_to_GCA_003402215.1.bigChain.link.bb?rlkey=08der7x6ap2pew89lb09gnn0v&dl=0) |
| `GCA_900093535.1_to_GCA_900093545.1.bigChain.bb` (→ PvT01) | [bb](https://www.dropbox.com/scl/fi/3jr31bup986nyn97al55t/GCA_900093535.1_to_GCA_900093545.1.bigChain.bb?rlkey=4c8woqts0y0uay4a0x31tun7k&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/9vd4shhdycgcsqc8znj28/GCA_900093535.1_to_GCA_900093545.1.bigChain.link.bb?rlkey=s3hwpkbs7ewokv4lrg20kgqiq&dl=0) |
| `GCA_900093535.1_to_GCA_040114635.1.bigChain.bb` (→ MHC087) | [bb](https://www.dropbox.com/scl/fi/6mwddcegzdq6yzyerxi4j/GCA_900093535.1_to_GCA_040114635.1.bigChain.bb?rlkey=j82jltmcrmrio9h33bc85kjij&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/5sss4ip8xceb0fk6ksvbm/GCA_900093535.1_to_GCA_040114635.1.bigChain.link.bb?rlkey=1jxmalwoubmy7zljgbprhksx5&dl=0) |

### MHC087 — `GCA_040114635.1/`

| File | Link |
|---|---|
| `trackDb.txt` | [Dropbox](https://www.dropbox.com/scl/fi/sl3ut6n4qy6zak1rzcsb7/trackDb.txt?rlkey=8jblxuwfn25s5bv01cw07p2my&dl=0) |
| `GCA_040114635.1.2bit` | [Dropbox](https://www.dropbox.com/scl/fi/riqza55y9lyod7nsas8ch/GCA_040114635.1.2bit?rlkey=g3fkkr42k5fl05f2apybvh6uv&dl=0) |
| `GCA_040114635.1_to_GCA_900093555.2.bigChain.bb` (→ PvP01) | [bb](https://www.dropbox.com/scl/fi/lxqmmsmkwrfk0ei889232/GCA_040114635.1_to_GCA_900093555.2.bigChain.bb?rlkey=6om9sg9g5wy5y3sdgaz6sxvu2&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/ney1tswpn7xtekfuivz3b/GCA_040114635.1_to_GCA_900093555.2.bigChain.link.bb?rlkey=utm8zkdv9w7fpkxdt0jcntnsz&dl=0) |
| `GCA_040114635.1_to_GCA_000002415.2.bigChain.bb` (→ Sal-I) | [bb](https://www.dropbox.com/scl/fi/18sqhcaxqyv8g5p5alb6h/GCA_040114635.1_to_GCA_000002415.2.bigChain.bb?rlkey=hgfq32j8h9kqae0wxd0zw2u0u&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/0u75v9zeapq29t6zlqy4k/GCA_040114635.1_to_GCA_000002415.2.bigChain.link.bb?rlkey=aa81drv8uepvbff5m7n3hnu2w&dl=0) |
| `GCA_040114635.1_to_GCA_914969965.1.bigChain.bb` (→ PvW1) | [bb](https://www.dropbox.com/scl/fi/yu1bt43io7udndkocduqm/GCA_040114635.1_to_GCA_914969965.1.bigChain.bb?rlkey=klxnq0lghyrzbm7n2cs587zmq&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/tm17yc8fedvlqcfve0uw8/GCA_040114635.1_to_GCA_914969965.1.bigChain.link.bb?rlkey=zd0t9085fx2wwrm7du81o0kmb&dl=0) |
| `GCA_040114635.1_to_GCA_949152365.1.bigChain.bb` (→ PAM) | [bb](https://www.dropbox.com/scl/fi/kun9x9u2oah5ksx1b6o9y/GCA_040114635.1_to_GCA_949152365.1.bigChain.bb?rlkey=aisdnjbskq6netzjomlcgpnl1&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/eg76xgzl9ygk6sexncjsx/GCA_040114635.1_to_GCA_949152365.1.bigChain.link.bb?rlkey=kwgtijhr1j0e6xy80drjx6nko&dl=0) |
| `GCA_040114635.1_to_GCA_003402215.1.bigChain.bb` (→ PvSY56) | [bb](https://www.dropbox.com/scl/fi/07w3u33kdj1rbl5fgvuyg/GCA_040114635.1_to_GCA_003402215.1.bigChain.bb?rlkey=4rcnp77zvaop2fweqtcrb61td&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/tixl7mrvdraxd9ia6qvsn/GCA_040114635.1_to_GCA_003402215.1.bigChain.link.bb?rlkey=t4tl2a8jpeyy75p5j0uhz2njb&dl=0) |
| `GCA_040114635.1_to_GCA_900093545.1.bigChain.bb` (→ PvT01) | [bb](https://www.dropbox.com/scl/fi/z82y07soil4qtcyinltxm/GCA_040114635.1_to_GCA_900093545.1.bigChain.bb?rlkey=tt0flgt6c9xdlribwqgk86h5f&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/ji5auueqiwwrbtjbukqh5/GCA_040114635.1_to_GCA_900093545.1.bigChain.link.bb?rlkey=b5ntmw8xfo99p06eovope999i&dl=0) |
| `GCA_040114635.1_to_GCA_900093535.1.bigChain.bb` (→ PvC01) | [bb](https://www.dropbox.com/scl/fi/gefgd6cvyib7fxjzt170c/GCA_040114635.1_to_GCA_900093535.1.bigChain.bb?rlkey=vw4e5rsje8oy3mr9km7xbgwsi&dl=0) · [link.bb](https://www.dropbox.com/scl/fi/81jvtm5ksndfob2jmefwx/GCA_040114635.1_to_GCA_900093535.1.bigChain.link.bb?rlkey=9lcircru8in22ciimehgpl8ym&dl=0) |

## Build scripts (in Git)

- [`tools/chain_to_bigChain.py`](https://github.com/nekrut/Pv4-pangenome/blob/main/v3/tools/chain_to_bigChain.py) — chain → bigChain.bed + bigLink.bed in one pass
- [`tools/bigChain.as`](https://github.com/nekrut/Pv4-pangenome/blob/main/v3/tools/bigChain.as) · [`tools/bigLink.as`](https://github.com/nekrut/Pv4-pangenome/blob/main/v3/tools/bigLink.as) — UCSC schema files for `bedToBigBed -as=…`

Run-time for the full 56-pair conversion: ~5 min on one CPU.
