const P_BITS: u64 = 57;
const MIN_LOOP: usize = 8;
const P57: u64 = 144115188075855859;
const TINYMT_MASK: u64 = 0x7fffffffffffffff;
const BYTES_UNDER_P: usize = 7;
const CHUNK_ALIGN: usize = 56;

fn main() {
    let (client_config, server_config) = init("abcdefghijklmnopqrstuvwxyz".as_bytes().to_vec());
}

struct ClientConfig {
    rows: usize,
    cols: usize,
    secret_m_vector: Vec<u64>,
    secret_n_vector: Vec<u64>,
}

struct ServerConfig {
    rows: usize,
    cols: usize,
    file: Vec<u8>,
}

struct tinymt64 {
    status: [u64; 2],
    mat1: u32,
    mat2: u32,
    tmat: u64,
}

impl tinymt64 {
    fn init(&mut self, seed: u64) {
        self.status[0] = seed ^ ((self.mat1 as u64) << 32);
        self.status[1] = self.mat2 as u64 ^ self.tmat;
        for i in 1..MIN_LOOP {
            self.status[i & 1] ^= (i as u128
                + 6364136223846793005_u128
                    * (
                        (self.status[(i - 1) & 1] ^ (self.status[(i - 1) & 1] >> 62)) as u128
                    )
                ) as u64;
        }
        self.period_certification();
    }

    fn period_certification(&mut self) {
        if (self.status[0] & TINYMT_MASK) == 0 && self.status[1] == 0 {
            self.status[0] = 'T' as u64;
            self.status[1] = 'M' as u64;
        }
    }

    fn rand_mod_p(&mut self) -> u64 {
        let mask = (1_u64 << P_BITS) - 1;
        let mut val: u64;
        loop {
            val = self.generate_uint64() & mask;
            if val < P57 {
                break;
            }
        }
        val
    }

    fn generate_uint64(&mut self) -> u64 {
        self.next_state();
        self.temper()
    }

    fn temper(&mut self) -> u64 {
        let mut x: u64;
        x = ((self.status[0] as u128 + self.status[1] as u128) % u64::MAX as u128) as u64;
        x ^= self.status[0] >> 8;
        x ^= (-((x & 1) as i64) & self.tmat as i64) as u64;
        x
    }

    fn next_state(&mut self) {
        let mut x: u64;
        self.status[0] &= TINYMT_MASK;
        x = self.status[0] ^ self.status[1];
        x ^= x << 12;
        x ^= x >> 32;
        x ^= x << 32;
        x ^= x << 11;
        self.status[0] = self.status[1];
        self.status[1] = x;
        self.status[0] ^= (-((x & 1) as i64) & self.mat1 as i64) as u64;
        self.status[1] ^= (-((x & 1) as i64) & ((self.mat2 as u64) << 32) as i64) as u64;
    }

    fn rand_vector(&mut self, size: usize) -> Vec<u64> {
        let mut vector = vec![0; size];
        for i in 0..size {
            vector[i] = self.rand_mod_p();
        }
        vector
    }
}

fn init(mut file: Vec<u8>) -> (ClientConfig, ServerConfig) {
    let num_chunks = 1 + (file.len() - 1) / BYTES_UNDER_P;
    let n =
        (((num_chunks as f64).sqrt() / CHUNK_ALIGN as f64).ceil() * CHUNK_ALIGN as f64) as usize;
    let m = 1 + (num_chunks - 1) / n;

    let seed = 2020;
    let mut state = tinymt64 {
        status: [0; 2],
        mat1: 0,
        mat2: 0,
        tmat: 0,
    };
    state.init(seed);

    let mut vector_u = vec![0; m];
    for i in 0..m {
        vector_u[i] = state.rand_mod_p();
    }

    let mut partials1 = vec![0_u128; n];
    let bytes_per_row = BYTES_UNDER_P * n;
    let chunk_mask = (1_u64 << (8 * BYTES_UNDER_P)) - 1;
    let file_extended: Vec<u8> = file.clone().into_iter().chain(vec![0; bytes_per_row]).collect();
    // file_extended.append(vec![0; bytes_per_row].as_mut());

    for i in 0..m {
        let mut raw_ind = 0;
        let raw_row = file_extended[(bytes_per_row * i)..bytes_per_row * (i + 1)].to_vec();
        let raw_row = raw_row.chunks_exact(8).map(|x| u64::from_le_bytes(x.try_into().unwrap())).collect::<Vec<u64>>();
        for full_ind in (0..n).step_by(8) {
            let mut data_val = (raw_row[raw_ind] & chunk_mask) as u128;
            partials1[full_ind] += data_val * vector_u[i] as u128;

            for k in 1..7 {
                let data_val = ((raw_row[raw_ind + k - 1] >> (64 - k * 8))
                    | ((raw_row[raw_ind + k] << (k * 8)) & chunk_mask)) as u128;
                partials1[full_ind + k] += data_val * vector_u[i] as u128;
            }

            data_val = (raw_row[raw_ind + 6] >> 8) as u128;
            partials1[full_ind + 7] += data_val * vector_u[i] as u128;

            raw_ind += 7;
        }
    }
    for k in 0..n {
        partials1[k] %= P57 as u128;
    }

    let client_config = ClientConfig{
        rows: n,
        cols: m,
        secret_m_vector: vector_u,
        secret_n_vector: partials1.iter().map(|x| *x as u64).collect::<Vec<u64>>(),
    };

    let server_config = ServerConfig{
        rows: n,
        cols: m,
        file: file,
    };

    return (client_config, server_config);
}

struct Client {
    config: ClientConfig,
}

impl Client {
    fn new(config: ClientConfig) -> Self {
        Self {
            config,
        }
    }

    // fn make_challenge_vector(&self, n: usize) -> Vec<u64> {
    //     let seed = 2020;

    // }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_init() {
        let (client_config, server_config) = init("abcdefghijklmnopqrstuvwxyz".as_bytes().to_vec());
        assert_eq!(client_config.rows, 56);
        assert_eq!(client_config.cols, 1);
        assert_eq!(client_config.secret_m_vector, vec![57829946736570845]);
        assert_eq!(client_config.secret_n_vector, vec![120891374367124132, 131035456404565768, 141179538442007404, 22254200296736808, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(server_config.rows, 56);
        assert_eq!(server_config.cols, 1);
        assert_eq!(server_config.file, "abcdefghijklmnopqrstuvwxyz".as_bytes().to_vec());
    }
}

// 			/* Audit */
// 			// send op code to server
// 			op = 'A';
// 			my_fwrite(&op, 1, 1, sock);
// 			fflush(sock);
			
// 			// create and send challenge vectors (of size n)
// 			start_time(&timer);				/* START COMP TIMER */
//             start_cpu_time(&cpu_timer);
// 			uint64_t* challenge1;
// 			uint64_t challengeBytes = n * sizeof(uint64_t);
// 			challenge1 = makeChallengeVector(n);
// 			client_comp_time = stop_time(&timer);		/* PAUSE COMP TIMER */
//             client_cpu_time = stop_cpu_time(&cpu_timer);
// 			start_time(&timer);				/* START COMM TIMER */
// 			my_fwrite(challenge1, 1, challengeBytes, sock);
// 			fflush(sock);

// 			// wait for ACK from server
// 			char ack = '0';
// 			my_fread(&ack, 1, 1, sock);
// 			if (ack == '1') comm_time = stop_time(&timer);	/* STOP COMM TIMER */
// 			else printf("Did not receive ACK from server after sending challenge.\n");
// 			printf("challenge[0] = %"PRIu64"\n", challenge1[0]);
// 			printf("challenge[n-1] = %"PRIu64"\n", challenge1[n-1]);

// 			// read response vectors from server (of size m)
// 			uint64_t* response1 = calloc(m, sizeof(uint64_t));
// 			uint64_t responseBytes = m * sizeof(uint64_t);
// 			my_fread(response1, 1, responseBytes, sock);
// 			printf("response[0] = %"PRIu64"\n", response1[0]);
// 			printf("response[m-1] = %"PRIu64"\n", response1[m-1]);

// 			// send previous comm_time as ack to server
// 			my_fwrite(&comm_time, sizeof(comm_time), 1, sock);
// 			fflush(sock);
// 			fprintf(stderr, "Sent 1-way comm time of %f to server.\n", comm_time);

// 			// run audit and report to client
// 			// use m for size
// 			start_time(&timer);				/* RESUME COMP TIMER */
//             start_cpu_time(&cpu_timer);
// 			int audit = runAudit(fconfig, challenge1,
// 							response1, n, m);
// 			client_comp_time += stop_time(&timer);		/* STOP TIMER */
//             client_cpu_time += stop_cpu_time(&cpu_timer);
// 			printf("Audit has ");
// 			printf(audit ? "PASSED!\n" : "FAILED.\n");

// 			//report computation time
// 			fprintf(stderr, "***CLIENT COMP TIME: %f***\n***CLIENT CPU  TIME: %f ***\n***CLIENT COMM TIME: %f ***\n", client_comp_time, client_cpu_time, comm_time);

// 			// clean up
// 			free(challenge1);
// 			free(response1);
// 			break;

//             uint64_t* makeChallengeVector(uint64_t size) {
// 	// seed Tiny Mersenne Twister
// 	uint64_t seed;
// #if __APPLE__
// 	if (getentropy(&seed, sizeof seed) == -1) {
// 		fprintf(stderr, "getentropy failed\n");
// 		exit(6);
// 	}
// #else
// 	if (getrandom(&seed, sizeof seed, 0) != sizeof(seed)) {
// 		fprintf(stderr, "getrandom failed\n");
// 		exit(7);
// 	}
// #endif
// 	tinymt64_t state = {0};
// 	tinymt64_init(&state, seed);

// 	// construct the randomized vector
// 	uint64_t* vector = calloc(size, sizeof(uint64_t));
// 	for (int i = 0; i < size; i++) {
// 		vector[i] = rand_mod_p(&state);
// 	}

// 	return vector;
// }

// int runAudit(FILE* fconfig, uint64_t* challenge1,
// 				uint64_t* response1, uint64_t n, uint64_t m) {
// 	uint128_t rxr1 = 0;
// 	uint128_t sxc1 = 0;

// 	// compute dot products:
// 	// random dot response & secret dot challenge.
// 	// config file read through once
// 	// doing modulo calc for each mul,
// 	// doing one modulo after all addition at end.
// 	size_t accum_count = 0;

// 	uint64_t *temp = malloc(MAX(m,n) * sizeof *temp);
// 	my_fread(temp, sizeof *temp, m, fconfig);

// 	for (int i = 0; i < m; i++) { /*random1 dot response1 (m)*/
// 		if ((accum_count += 2) > MAX_ACCUM_P) {
// 			rxr1 %= P57;
// 			accum_count = 2;
// 		}
// 		rxr1 += ((uint128_t)temp[i]) * response1[i];
// 	}
// 	rxr1 %= P57;
// 	accum_count = 0;
// 	my_fread(temp, sizeof *temp, n, fconfig);
// 	for (int i = 0; i < n; i++) { /*secret1 dot challenge1 (n)*/
// 		if ((accum_count += 2) > MAX_ACCUM_P) {
// 			sxc1 %= P57;
// 			accum_count = 2;
// 		}
// 		sxc1 += ((uint128_t)temp[i]) * challenge1[i];
// 	}
// 	sxc1 %= P57;
// 	free(temp);

// 	// check for equal and return result
// 	// 1 for pass
// 	// 0 for fail (default)
// 	printf("rxr1 = %"PRIu64"\n", (uint64_t)rxr1);
// 	printf("sxc1 = %"PRIu64"\n", (uint64_t)sxc1);
// 	return (rxr1 == sxc1);
// }
