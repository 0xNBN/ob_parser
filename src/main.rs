use std::time::Instant;
use std::hint::black_box;
use std::arch::x86_64::*;
use std::ptr;

fn main() {

    if cfg!(target_arch = "x86_64") {
        println!("CPU supports:");

        if std::is_x86_feature_detected!("sse") {
            println!("  - SSE");
        }
        if std::is_x86_feature_detected!("sse2") {
            println!("  - SSE2");
        }
        if std::is_x86_feature_detected!("avx") {
            println!("  - AVX");
        }
        if std::is_x86_feature_detected!("avx2") {
            println!("  - AVX2");
        }
        if std::is_x86_feature_detected!("avx512f") {
            println!("  - AVX-512");
        }
    } else {
        println!("Not running on an x86_64 architecture.");
    }


    //let message = br#"{"stream":"btcusdt@bookTicker","data":{"u":55361077872,"s":"BTCUSDT","b":"94906.44000000","B":"2.61478000","a":"94906.45000000","A":"2.43671000"}}"#;
    let message = br#"{"stream":"btcusdt@bookTicker","data":{"u":55361077872,"s":"BTCUSDT","b":"94906.44000000","B":"2.61478000","a":"94906.45000000","A":"2.43671000"}}"#;
    
    let mut times = vec![];
    for _ in 0..100000 {
        let now = Instant::now();
        let b = parse_bookticker(message);
        times.push(now.elapsed().as_nanos() as u64);
        //println!("{:?}",b);
        black_box(b);
    }
    
    let s:u64 = times[50000..].iter().sum();
    println!("{:?}", s / 50000 as u64);
}



fn parse_bookticker(input: &[u8]) -> Result<(usize, (u64, u64), (u64, u64)), Box<dyn std::error::Error>> {
   

    let ptr = &input.as_ptr(); 
    let len_tot=input.len();

    // Searching for the indexes of the things to parse
    
    let start_u=43;
    let mut end_u=0;
    //Starting at 54 because we assume the ob update id will only go increasing
    // Meaning we'll always have at least 11 digits
    for i in 54..len_tot{

        if input[i]==44{
            end_u=i;
            break
        }
    }

    // Discarding static data : pair name
    // We might want to have the static sequence as an input to handle multiple pairs
    let start_b1=end_u+20;
    let mut end_b1=0;
    for i in start_b1..len_tot{
        if input[i]==46{
            end_b1=i;
            break
        }
    }
    let start_b2=end_b1+1;
    let end_b2=start_b2+8;

    let start_B1=end_b2+7;
    let mut end_B1=0;
    for i in start_B1..len_tot{
        if input[i]==46{
            end_B1=i;
            break
        }
    }
    let start_B2=end_B1+1;
    let end_B2=start_B2+8;

    let start_a1=end_B2+7;
    let mut end_a1=0;
    for i in start_a1..len_tot-1{
        
        if input[i]==46{
            end_a1=i;
            break
        }
    }

    // For the last field we avoid one search because the end is deterministic : 
    // because } and fixed number of decimals so we know whats the index of the last "."
    let start_a2=end_a1+1;
    let end_a2=start_a2+8;

    let start_A1=end_a2+7;
    let end_A2=len_tot-3;
    let start_A2= end_A2-8;
    let end_A1=start_A2-1;


    // Lets start the parsing
    let mut uu:u64=0;
    let mut bb:u64=0;
    // This is the decimal pard of bb2
    let mut bb2:u32=0;
    let mut BB:u64=0;
    //Decimal part of BB
    let mut BB2:u32=0;
    let mut aa:u64=0;
    // Decimal part
    let mut aa2:u32=0;
    let mut AA:u64=0;
    // Decimal part
    let mut AA2:u32=0;


    unsafe {


        // PARSING U : we assume it'll be less than 16 digits

        //11 bytes here max 16
        let part1 = &input[start_u..end_u];
        let len_u =end_u-start_u;
        // Load 128 bits (16 bytes) starting from the pointer into an SIMD register
        
        let data = _mm_loadu_si128(part1.as_ptr() as *const __m128i);

        // Subtract '0' from each byte to convert ASCII characters to digits
        let zero = _mm_set1_epi8(b'0' as i8);
        let digits = _mm_sub_epi8(data, zero);
        // Process only the required bytes (e.g., first 11 digits)
        let mut result = [0u8; 16];
        _mm_storeu_si128(result.as_mut_ptr() as *mut __m128i, digits);

        for i in 0..len_u{
            uu = (uu << 3) + (uu << 1) + (result[i] as u64);
        }
        














        //----------------PARSING THE 2 DECIMAL PARTS OF b and a in 1 SIMD register----------


        // We know it will be 16 digits as decimals as 8 digits for all numbers
        let part1 = &input[start_b2..end_b2];
        let part2 = &input[start_a2..end_a2];


        // Create a single contiguous byte array (pre-allocate if necessary)
        let mut combined: [u8; 16] = [0; 16];

        // Manually copy the bytes into the array
        ptr::copy_nonoverlapping(part1.as_ptr(), combined.as_mut_ptr(), 8);
        ptr::copy_nonoverlapping(part2.as_ptr(), combined.as_mut_ptr().add(8), 8);
        let data = _mm_loadu_si128(combined.as_ptr() as *const __m128i);

        // Subtract '0' from each byte to convert ASCII characters to digits
    
        let digits = _mm_sub_epi8(data, zero);
        let mut result = [0u8; 16];
        _mm_storeu_si128(result.as_mut_ptr() as *mut __m128i, digits);

  
        // Fixed 8 digits so we unroll

        bb2 = ((bb2 << 3) + (bb2 << 1)) + (result[0]  as u32);
        bb2 = ((bb2 << 3) + (bb2 << 1)) + (result[1]  as u32);
        bb2 = ((bb2 << 3) + (bb2 << 1)) + (result[2]  as u32);
        bb2 = ((bb2 << 3) + (bb2 << 1)) + (result[3]  as u32);
        bb2 = ((bb2 << 3) + (bb2 << 1)) + (result[4]  as u32);
        bb2 = ((bb2 << 3) + (bb2 << 1)) + (result[5]  as u32);
        bb2 = ((bb2 << 3) + (bb2 << 1)) + (result[6]  as u32);
        bb2 = ((bb2 << 3) + (bb2 << 1)) + (result[7]  as u32);

      

        aa2 = ((aa2 << 3) + (aa2 << 1)) + (result[8] as u32);
        aa2 = ((aa2 << 3) + (aa2 << 1)) + (result[9] as u32);
        aa2 = ((aa2 << 3) + (aa2 << 1)) + (result[10] as u32);
        aa2 = ((aa2 << 3) + (aa2 << 1)) + (result[11] as u32);
        aa2 = ((aa2 << 3) + (aa2 << 1)) + (result[12] as u32);
        aa2 = ((aa2 << 3) + (aa2 << 1)) + (result[13] as u32);
        aa2 = ((aa2 << 3) + (aa2 << 1)) + (result[14] as u32); 
        aa2 = ((aa2 << 3) + (aa2 << 1)) + (result[15] as u32); 









        //-------------PARSING THE 2 DECIMAL PARTS OF B and A in 1 SIMD register-----------




        let part1 = &input[start_B2..end_B2];
        let part2 = &input[start_A2..end_A2];


        // Create a single contiguous byte array (pre-allocate if necessary)
        let mut combined: [u8; 16] = [0; 16];

        // Manually copy the bytes into the array
        ptr::copy_nonoverlapping(part1.as_ptr(), combined.as_mut_ptr(), 8);
        ptr::copy_nonoverlapping(part2.as_ptr(), combined.as_mut_ptr().add(8), 8);
        let data = _mm_loadu_si128(combined.as_ptr() as *const __m128i);

        // Subtract '0' from each byte to convert ASCII characters to digits
    
        let digits = _mm_sub_epi8(data, zero);
        let mut result = [0u8; 16];
        _mm_storeu_si128(result.as_mut_ptr() as *mut __m128i, digits);



        BB2 = ((BB2 << 3) + (BB2 << 1)) + (result[0]  as u32);
        BB2 = ((BB2 << 3) + (BB2 << 1)) + (result[1]  as u32);
        BB2 = ((BB2 << 3) + (BB2 << 1)) + (result[2]  as u32);
        BB2 = ((BB2 << 3) + (BB2 << 1)) + (result[3]  as u32);
        BB2 = ((BB2 << 3) + (BB2 << 1)) + (result[4]  as u32);
        BB2 = ((BB2 << 3) + (BB2 << 1)) + (result[5]  as u32);
        BB2 = ((BB2 << 3) + (BB2 << 1)) + (result[6]  as u32);
        BB2 = ((BB2 << 3) + (BB2 << 1)) + (result[7]  as u32);

    

        AA2 = ((AA2 << 3) + (AA2 << 1)) + (result[8] as u32);
        AA2 = ((AA2 << 3) + (AA2 << 1)) + (result[9] as u32);
        AA2 = ((AA2 << 3) + (AA2 << 1)) + (result[10] as u32);
        AA2 = ((AA2 << 3) + (AA2 << 1)) + (result[11] as u32);
        AA2 = ((AA2 << 3) + (AA2 << 1)) + (result[12] as u32);
        AA2 = ((AA2 << 3) + (AA2 << 1)) + (result[13] as u32);
        AA2 = ((AA2 << 3) + (AA2 << 1)) + (result[14] as u32); 
        AA2 = ((AA2 << 3) + (AA2 << 1)) + (result[15] as u32); 



















        //--------- PARSING THE AMOUNTS WITHOUT DECIMAL FOR b and a
        // If we can encapsulate the 2 amounts into 1 same SIMD we do it
        // Otherwise we parse the 2 amounts independantly inside a SIMD
        // This works for amounts of tokens up to 1 million billions (16 digits without the decimals)
        if end_a1-start_a1+end_b1-start_b1<=16{

            let part1 = &input[start_b1..end_b1];
            let part2 = &input[start_a1..end_a1];
            let len_a1=end_a1-start_a1;
            let len_b1=end_b1-start_b1;

            // Create a single contiguous byte array (pre-allocate if necessary)
            let mut combined: [u8; 16] = [0; 16];

            // Manually copy the bytes into the array
            ptr::copy_nonoverlapping(part1.as_ptr(), combined.as_mut_ptr(), len_b1);
            ptr::copy_nonoverlapping(part2.as_ptr(), combined.as_mut_ptr().add(len_b1), len_a1);
            let data = _mm_loadu_si128(combined.as_ptr() as *const __m128i);

            // Subtract '0' from each byte to convert ASCII characters to digits
        
            let digits = _mm_sub_epi8(data, zero);
            let mut result = [0u8; 16];
            _mm_storeu_si128(result.as_mut_ptr() as *mut __m128i, digits);

            // Process only the required bytes (e.g., first 11 digits)


            for i in 0..len_b1{
                bb = ((bb << 3) + (bb << 1)) + (result[i]  as u64);
            }

            for i in len_b1..(len_b1+len_a1){
                aa = ((aa << 3) + (aa << 1)) + (result[i-len_b1]  as u64);
            }
       
            // Adding decimals -> Adapting non decimal part 

            for i in 0..8{
                bb = ((bb << 3) + (bb << 1));
                aa = ((aa << 3) + (aa << 1));
            }
        
            bb=bb+bb2 as u64;
            aa=aa+aa2 as u64;
     



        }// Else part to implement : 2 independant SIMD










        //--------- PARSING THE AMOUNTS WITHOUT DECIMAL FOR B and A
        // If we can encapsulate the 2 amounts into 1 same SIMD we do it
        // Otherwise we parse the 2 amounts independantly inside a SIMD
        // This works for amounts of tokens up to 1 million billions (16 digits without the decimals)
        if end_A1-start_A1+end_B1-start_B1<=16{

            let part1 = &input[start_B1..end_B1];
            let part2 = &input[start_A1..end_A1];
            let len_A1=end_A1-start_A1;
            let len_B1=end_B1-start_B1;


            // Create a single contiguous byte array (pre-allocate if necessary)
            let mut combined: [u8; 16] = [0; 16];

            // Manually copy the bytes into the array
            ptr::copy_nonoverlapping(part1.as_ptr(), combined.as_mut_ptr(), len_B1);
            ptr::copy_nonoverlapping(part2.as_ptr(), combined.as_mut_ptr().add(len_B1), len_A1);
            let data = _mm_loadu_si128(combined.as_ptr() as *const __m128i);

            // Subtract '0' from each byte to convert ASCII characters to digits
        
            let digits = _mm_sub_epi8(data, zero);
            let mut result = [0u8; 16];
            _mm_storeu_si128(result.as_mut_ptr() as *mut __m128i, digits);

            // Process only the required bytes (e.g., first 11 digits)


            for i in 0..len_B1{
                BB = ((BB << 3) + (BB << 1)) + (result[i]  as u64);
            }

            for i in len_B1..(len_B1+len_A1){
                AA = ((AA << 3) + (AA << 1)) + (result[i-len_B1]  as u64);
            }
           
            // Adding decimals -> Adapting non decimal part 
            for i in 0..8{
                BB = ((BB << 3) + (BB << 1));
                AA = ((AA << 3) + (AA << 1));
            }
           
            BB=BB+BB2 as u64;
            AA=AA+AA2 as u64;
            
        }





    }
    



    // Ensure all fields are parsed

    Ok((uu as usize, (bb as u64, BB), (aa as u64, AA)))

}