window.BENCHMARK_DATA = {
  "lastUpdate": 1725984055365,
  "repoUrl": "https://github.com/manozzing/EnvironmentalTransport.jl",
  "entries": {
    "Julia benchmark result": [
      {
        "commit": {
          "author": {
            "email": "ctessum@gmail.com",
            "name": "Christopher Tessum",
            "username": "ctessum"
          },
          "committer": {
            "email": "ctessum@gmail.com",
            "name": "Christopher Tessum",
            "username": "ctessum"
          },
          "distinct": true,
          "id": "a93ea925a4ddab2621d63a6d6f69c97764dac8d5",
          "message": "Fix typo",
          "timestamp": "2024-09-04T13:22:31-05:00",
          "tree_id": "8c24152791605434d5f1f5fe2a68c49a08c56da7",
          "url": "https://github.com/manozzing/EnvironmentalTransport.jl/commit/a93ea925a4ddab2621d63a6d6f69c97764dac8d5"
        },
        "date": 1725984055018,
        "tool": "julia",
        "benches": [
          {
            "name": "Advection Simulator/in-place/l94_stencil/0.625 x 0.5 (N=31719)",
            "value": 89690660,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4416\nallocs=19\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Advection Simulator/in-place/l94_stencil/0.3125 x 0.25 (N=126222)",
            "value": 352593588,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4416\nallocs=19\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Advection Simulator/in-place/ppm_stencil/0.625 x 0.5 (N=31719)",
            "value": 9423476648,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4416\nallocs=19\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Advection Simulator/in-place/ppm_stencil/0.3125 x 0.25 (N=126222)",
            "value": 37418523124,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4416\nallocs=19\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Advection Simulator/out-of-place/l94_stencil/0.625 x 0.5 (N=31719)",
            "value": 91599505.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=257552\nallocs=21\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Advection Simulator/out-of-place/l94_stencil/0.3125 x 0.25 (N=126222)",
            "value": 352756177,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1013584\nallocs=21\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Advection Simulator/out-of-place/ppm_stencil/0.625 x 0.5 (N=31719)",
            "value": 9573658931,
            "unit": "ns",
            "extra": "gctime=0\nmemory=257552\nallocs=21\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Advection Simulator/out-of-place/ppm_stencil/0.3125 x 0.25 (N=126222)",
            "value": 37353787118,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1013584\nallocs=21\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      }
    ]
  }
}