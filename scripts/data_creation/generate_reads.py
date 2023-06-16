import random
import uuid
import string


READ_COUNT = 40

LEN_A_TAIL = 40
LEN_REMAINING_SEQ = 60

LEN_TOTAL = LEN_A_TAIL + LEN_REMAINING_SEQ
VALID_PHRED = [chr(i) for i in range(33,127)]

for i in range(40):
    print(f"@{uuid.uuid4()}")
    random_seq = "".join([random.choice("ACTG") for _ in range(LEN_REMAINING_SEQ)])
    print("A" *LEN_A_TAIL + random_seq)
    print("+")
    print("".join([random.choice(VALID_PHRED) for _ in range(LEN_TOTAL)]))