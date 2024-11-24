# symmetric-key-threshold

## Setup
See `MP-SPDZ/README.md` for instructions on setting up the MP-SPDZ library.

## Running MPC 
```
cd MP-SPDZ
```
Compilation
```
/compile.py -P 57896044618658097711785492504343953926634992332820282019728792003956566065153 mimc [n_clients] && PLAYERS=[n_clients] Scripts/mascot.sh mimc-[n_clients] --batch-size 20
```

In separate terminal
```
python3 main.py
```

To run the verifier 
- Run `main.py` but set `TESTING=True` (This writes outputs to files)
```
python3 verifier.py
```
