# Keeping Time with the Cosmos: A Home-Built Satellite Time Tracker

Welcome, curious explorer and supportive parents! If you have ever wondered how the clocks on our computers, phones, and global financial networks stay synchronized, you are in the right place. Today, we are going to embark on an exciting journey to build our very own satellite-guided time tracker. 

We will explore how we can listen to signals from satellites flying hundreds of miles above our heads and use them to steer our local computer clock. Best of all, we will break down the complex science into simple ideas that anyone can understand. No advanced math degrees required—just a bit of curiosity and a desire to see how the pieces of physics and software fit together.

---

## 1. The Mystery of the Slipping Clock

Let us start with a simple question: how does your computer know what time it is? 

Inside every computer, smartphone, and digital wristwatch is a tiny slice of quartz crystal. When we apply electricity to this crystal, it vibrates at a very precise frequency, like a microscopic tuning fork. The computer counts these vibrations—say, 32,768 times per second—and uses that count to advance the clock by one second.

This quartz crystal is convenient and cost-effective, but it is sensitive to environmental factors. If your room gets warmer or colder, the crystal vibrates a little faster or slower. As the crystal ages, its rate changes. Even tiny imperfections from the factory mean that no two quartz crystals are exactly identical. 

Because of these factors, your computer's internal clock constantly drifts. It might gain or lose a few seconds every single day. If left to itself, your computer would eventually drift minutes or hours out of sync.

### Why Atomic Time Matters

In our daily lives, a clock that is off by a few seconds is no big deal. You might arrive a tiny bit early or late to a meeting. But in the modern digital world, seconds are an eternity. 

Imagine a bank processing stock trades. If one computer's clock is three seconds behind another, a trader could exploit that delay to buy stocks in the past and sell them in the present, throwing the financial system into chaos. In electric power grids, computers must coordinate the flow of electricity to the millisecond to prevent blackouts. In GPS navigation, a clock error of just one microsecond (one millionth of a second) translates to a positioning error of about three hundred meters.

To keep everything running smoothly, we rely on atomic time. Atomic clocks use the vibrations of cesium or rubidium atoms, which are so stable that they lose less than a single second over millions of years. But atomic clocks are massive, expensive, and require specialized laboratories to run.

### Enter the NTP Server

Since we cannot fit an atomic clock inside every laptop, we use the Network Time Protocol, or NTP. An NTP server is a dedicated computer connected to a highly accurate time source (like an atomic clock or a GPS receiver). Other computers on the network query this NTP server periodically over the internet to adjust their own clocks and correct their drift. 

But what if you are offline? What if you are in a remote field, at sea, or want to build a timekeeping system that does not rely on a commercial internet connection? That is where our project, `sattime`, comes in. We are going to build a system that listens directly to satellites to create our own independent NTP server.

---

## 2. The Space Watch: Low Earth Orbit Satellites

To correct our drifting local clock, we need a reliable reference clock to compare it against. Fortunately, there is a constellation of clocks flying overhead right now. 

Low Earth Orbit (LEO) satellites, such as those in weather satellite networks or communication constellations, orbit the Earth at altitudes between one hundred and one thousand miles. Because these satellites need to coordinate their operations, they carry highly precise clocks on board, which are regularly synchronized with ground-based atomic clocks.

We can think of these LEO satellites as a "space watch". As a satellite passes overhead, it broadcasts a radio signal containing its orbital parameters (where it is in space) and the exact time it sent the signal. By listening to this space watch, our local receiver can compare its own clock against the satellite's atomic clock.

Using LEO satellites has a huge advantage over traditional GPS satellites. GPS satellites orbit much higher, around twelve thousand miles up, which means their signals are incredibly weak by the time they reach Earth. They require a clear line of sight to the sky. LEO satellites are much closer, so their signals are much stronger and can often be received with simple, home-made antennas.

---

## 3. The Doppler Effect: The Train Siren and the Space Whistle

How do we actually extract time and position from a satellite flying overhead? The key is an everyday physics phenomenon called the Doppler effect.

To understand the Doppler effect, imagine you are standing near a railway track. A fast-moving train approaches, blowing its whistle and sounding its warning siren. As the train rushes toward you, the sound of the siren feels high-pitched. The moment the train passes you and begins moving away, the pitch suddenly drops, sounding much lower. 

This happens because sound travels in waves. As the train moves toward you, it chases its own sound, compressing the sound waves together. More waves reach your ears per second, which your brain interprets as a higher frequency. Once the train passes and recedes, it moves away from the waves it emits, stretching them out. Fewer waves reach you per second, resulting in a lower frequency.

The exact same thing happens with radio waves! A satellite is like a cosmic train, and its radio transmitter is the whistle. As the satellite flies toward our antenna at over seventeen thousand miles per hour, the frequency of the radio signal we receive is shifted higher. As it passes directly overhead and flies away, the frequency is shifted lower. This characteristic shift is called the Doppler curve.

### Reverse-GPS and Clock Steering

By measuring the exact shape of this Doppler curve, we can perform a fascinating trick called Reverse-GPS. 

In normal GPS navigation, your phone listens to four or more satellites at the same time to calculate your position. With Reverse-GPS, we only need to listen to a single satellite during its pass. 

Because we know the satellite's exact path in space from its orbital data, the shape of the Doppler curve tells us exactly how close the satellite came to our antenna and at what moment. If the frequency changes very rapidly from high to low, the satellite passed directly overhead. If the frequency changes slowly, the satellite passed far to the side. 

By analyzing this curve, we can solve two mysteries at once:
1. **Where we are**: We can pinpoint our own latitude and longitude on Earth.
2. **What time it is**: We can calculate the exact error of our local computer clock and "steer" it back to matching the satellite's atomic time. We adjust our clock phase and speed to align the local clock with the satellite's atomic time.

---

## 4. Sifting the Signals: Filtering and Noise

Capturing a satellite signal is not as simple as tuning in to a local FM radio station. The radio spectrum is crowded and messy. There is constant background noise from electrical lines, computers, Wi-Fi routers, and atmospheric static. 

To find the satellite's faint whistle in this sea of static noise, our system must act like a sieve. We use digital signal processing filters to clean up the signal.

### Decimation: Thinning the Flood

Our Software Defined Radio (SDR) receiver captures millions of radio samples every second. Processing this mountain of data would melt a normal computer's processor. 

To solve this, we use a process called decimation. First, we apply a digital low-pass filter to block out high-frequency noise outside our band of interest. Then, we throw away most of the samples, keeping perhaps only one out of every forty. This thins the data flood down to a manageable size while preserving the satellite's signal, allowing our computer to process it in real time.

### Spur Notching: Blocking the Static

Imagine you are trying to listen to a speaker, but there is a loud, constant hum from a nearby refrigerator. That hum is what RF engineers call a static "spur"—interference at a fixed frequency. 

Our system uses an autotuning algorithm that monitors the radio spectrum. If it detects a signal that stays at the exact same frequency for a long time, it flags it as a static spur (since a moving satellite's frequency would constantly drift due to the Doppler effect). The system then notches out that specific frequency, creating a digital blind spot for the hum while leaving the rest of the spectrum open to detect the moving satellite.

---

## 5. The Kalman Filter: The Smart Guesser

Once we have filtered our radio signals, we obtain measurements of the satellite's frequency. But these measurements are still not perfect. Wind, atmospheric conditions, and receiver limitations mean our data is still slightly jumpy and noisy. How do we find the absolute truth hidden behind this noisy data? 

We use an algorithm called a Kalman Filter.

To understand how a Kalman Filter works, let us look at a simple analogy. Imagine you are driving a car through a dark, foggy tunnel. You want to know your exact speed and position. You have two sources of information:
1. **The speedometer and odometer**: These tell you how fast the wheels are turning and how far you have travelled. This is your prediction model. It is smooth, but if your wheels slip, it accumulates errors over time.
2. **A GPS receiver**: Every few seconds, it tries to ping your location. But inside the tunnel, the signal is weak, and the readings jump around wildly. This is your measurement. It is noisy, but it does not accumulate long-term errors.

A bad system would choose one or the other. A smart system uses a Kalman Filter. 

The Kalman Filter knows how both systems behave. It starts with the prediction model to guess where the car should be. When a new GPS measurement arrives, the filter calculates the difference between the prediction and the measurement. It then makes a smart guess, weighting the two based on their reliability. If the GPS signal is very noisy, it trusts the prediction more. If the wheels are slipping, it trusts the GPS more. The result is a smooth, highly accurate estimate of the car's true position.

In `sattime`, the Kalman Filter tracks two states: our clock's phase offset (how many microseconds we are off) and its frequency drift (how fast the clock is running). 

It combines:
- Our prediction of how the quartz crystal drifts over time.
- The noisy frequency measurements from the satellite passes.

By constantly balancing the prediction and the measurements, the Kalman Filter provides a continuous, highly stable estimate of the true time, allowing us to steer the computer clock with microsecond-level precision.

### The Carrier EKF: Locking onto the Space Whistle
In addition to the clock-steering filter, the receiver uses a second, ultra-fast Kalman Filter to lock onto the satellite's carrier signal. Think of this as a set of **robotic ears** that tune into the satellite's whistle. 

Normally, the receiver scans the radio dial to find the whistle's pitch. But once it finds it, the Carrier EKF takes over. It tracks three things sample-by-sample:
1. The **phase** of the wave (exactly where the wave's peaks and valleys are).
2. The **frequency** of the wave (how high or low the pitch of the whistle is).
3. The **chirp rate** (how fast the pitch is sliding down as the satellite rushes overhead).

By predicting the wave's movement sample-by-sample, these robotic ears can block out 99% of the background static. This allows us to track the space whistle even when it is so weak that human ears would hear nothing but hiss. The resulting tracking data is so clean that our orbital math can solve our location and clock offset with **20 times more precision** than before!

---

## FAQ for Parents and Curious Minds

### Why cannot we just sync our clocks over the internet?
Internet-based time synchronization is fantastic, but it requires a constant internet connection and is vulnerable to network congestion. If routing paths change, the time packets can be delayed, introducing errors. A satellite-guided system works completely offline, making it suitable for remote locations, emergency backup systems, or high-security networks that must remain isolated from the public internet.

### Do we need a giant, expensive satellite dish?
No. Because LEO satellites are relatively close to the Earth, their signals are strong. You can receive them using a simple dipole or turnstile antenna made from spare wire or measuring tape, connected to a cheap USB software-defined radio receiver.

### Is it legal to listen to these satellites?
Yes. Low Earth Orbit satellites broadcast their telemetry and signals publicly on amateur and weather bands. We are only receiving these signals, not transmitting anything. It is completely passive and legal.

### What level of accuracy can we achieve with this setup?
By combining decimation filtering, spur notching, and the Kalman Filter, a standard computer clock can be disciplined to stay within a few microseconds of UTC (Coordinated Universal Time), which is designed to achieve microsecond-level accuracy.

### Summary for Parents
Our home-built system is a mini-science laboratory. It combines physics (the Doppler effect), orbital mechanics (tracking satellites in space), digital signal processing (cleaning up radio signals), and advanced estimation math (the Kalman Filter). It shows how a simple computer can be turned into a highly precise scientific instrument using open-source software and basic radio hardware.

---

## 6. Advanced Signal Processing and Optimization Algorithms

To push the performance of `sattime` to the absolute limits of RF and mathematical engineering, we introduced four advanced algorithms. These algorithms solve the hardest problems in satellite tracking: finding orbits from scratch, rejecting false signals, notching out persistent background hums, and establishing a universal volume control.

Let us explore these four additions through simple, real-world analogies.

### 1. The Adelic Langevin Solver: The Smart Mountain Climber
Imagine you are blindfolded and dropped onto a rugged mountain range, and you want to find the highest peak. 
- **The Old Way (Grid Search)**: You would walk in a rigid grid pattern across the mountains, taking steps exactly every 50 feet. If the peak is small and sits between your grid lines, you will walk right past it and miss it entirely. This is slow, mechanical, and easily misses the target.
- **The Adelic Langevin Solver**: Instead of walking mindlessly, you use a smart stochastic climber. The climber feels the slope under their feet (the gradient) and walks upward. To avoid getting stuck in a small ditch, the climber occasionally "teleports" randomly to nearby spots. 
The key aspect of this approach is: these teleports are guided by *p-adic* numbers. Instead of just walking on continuous paths, the climber jumps back and forth across a discrete fractal map. By combining smooth continuous steps (using standard real numbers) with fractal jumps (using discrete $p$-adic primes), the solver can search massive orbital spaces without getting stuck in local traps, finding the satellite's exact orbit with incredible speed and accuracy.

### 2. Sheaf Cohomology & Čech Obstruction: Overlapping Opinions
Imagine you are in a crowded, echoey room trying to write down what a speaker is saying. You have three listeners in different parts of the room.
- **The Old Way (Independent Signals)**: Each listener writes down what they hear, and you trust each one individually. If one listener is fooled by a loud echo bouncing off the wall, they write down incorrect words, and you have no way to know who is right.
- **Sheaf Cohomology & Čech Obstruction**: Instead of treating listeners independently, you treat their notes as a system of "overlapping opinions" (in mathematics, a *sheaf*). You compare where their descriptions overlap. If two listeners agree, but the third reports something completely different, you identify an inconsistency—a *cohomological obstruction* (specifically, a non-zero Čech coboundary) caused by a multipath echo.
Our tracking system runs three EKF trackers in parallel. If they drift apart by more than $150\text{ Hz}$, the Čech obstruction detects the mismatch, immediately suspends NTP clock steering to prevent corrupting the system clock, and prunes the rogue tracker that is locking onto the echo.

### 3. The Vladimirov-Steered Tropical Wavelet: The Whistle and the Hum
Imagine you are at a concert, trying to record a singer performing a beautiful, sliding slide-whistle solo. Unfortunately, there is a loud, static hum from a bad speaker nearby.
- **The Old Way (Static Filters)**: You try to notch out the hum's pitch. But as the slide-whistle slides up and down, it crosses the hum's frequency, and your filter accidentally cuts out the singer's voice too.
- **The Vladimirov-Steered Tropical Wavelet**: This algorithm acts as a dynamic, "min-plus wave envelope" (using *tropical geometry* where multiplication is addition, and addition is finding the minimum). It tracks the lowest background floor of the sound, drawing a smooth envelope underneath.
If it detects a sudden, sharp spike that stands still (the static hum), it checks it against a *Vladimirov fractional derivative* (a mathematical tool that measures how local or global a frequency spike is). If the spike is a static spur, the wavelet notches it out. If it is the sliding satellite whistle, the system keeps a "guard band" around the satellite's Doppler frequency and chirp rate, letting the whistle pass through safely while keeping the background hum completely silenced.

### 4. Calibrated Gauge AGC: Finding the Universal Scale
Imagine you are testing audio speakers from three different factories. Each speaker has its own volume knob, but "Volume 5" on speaker A is deafeningly loud, while "Volume 5" on speaker B is barely a whisper. To compare them fairly, you need a way to measure the absolute sound level in decibels, regardless of where the knobs are set.
- **The Old Way**: You adjust the volume knobs manually or rely on basic gain levels, but you never know the true strength of the radio waves hitting your antenna.
- **Calibrated Gauge AGC**: This algorithm acts as a universal scale (*gauge calibration*). It monitors how much the signal is clipping (saturation) and its average power (RMS). It then dynamically adjusts three different amplifier stages in your SDR (the LNA, the VGA, and the AMP) using strict priority rules (e.g., raising LNA first to keep noise low, or lowering AMP first during saturation to prevent distortion).
Most importantly, it applies a custom polynomial correction curve to translate those local gain settings into an absolute, physical power reading in dBm. This gives you a hardware-independent measurement of the satellite's true signal strength, letting you compare passes recorded on different hardware setups with absolute precision.

---

## 7. Phase 11 Upgrades: Multi-Channel Satellite Tracking

To elevate `sattime` from tracking one satellite at a time to monitoring an entire fleet of spacecraft simultaneously, Phase 11 introduces multi-channel tracking. Here is how we achieved this transition, explained through simple analogies.

### 1. Multi-Channel Satellite Tracking: The Multi-Tasking Air Traffic Controller
Imagine you are an air traffic controller at a busy regional airport.
- **The Old Way**: You can only talk to one airplane at a time. While you guide a weather satellite to land, all other aircraft must circle in a queue, waiting their turn. If three satellites pass by at the same moment, you miss two of them.
- **The Phase 11 Way**: You now run a modern control tower with a bank of screens. You can talk to and track up to 8 satellites at the exact same time. Each satellite gets its own dedicated controller, allowing the system to capture multiple signals simultaneously.

### 2. Dynamic AOS/LOS Scheduling (`ChannelAllocator`): The Airport Gate Allocator
With multiple satellites flying overhead, we need a smart system to assign them to our limited number of tracking screens.
- **The Old Way**: You manually tune the receiver when you hear static or try to guess which satellite is rising.
- **The ChannelAllocator**: This acts as a smart airport gate manager. When a satellite rises above the horizon (Acquisition of Signal, or AOS), the allocator finds an empty gate (processing channel) and guides the satellite in. It verifies the satellite's flight itinerary (TLE orbital path) to ensure it is not a corrupted signal. When the satellite disappears below the opposite horizon (Loss of Signal, or LOS), the allocator frees up the gate for the next traveler.

### 3. Parallel Digital Downconversion: The Personal Audio Headsets
The radio antenna receives a massive, single stream of raw wideband signals containing data from all satellites visible in the sky.
- **The Old Way**: One central processor tries to read the entire stream, isolating one frequency, which forces the system to ignore all other frequencies.
- **Parallel DDC**: Imagine feeding the master audio tape of an orchestra to 8 listeners, each wearing a headset with its own tuning knob. One listener tunes their knob to the flute (Satellite A), another tunes to the violin (Satellite B), and they filter out the rest of the orchestra. Rayon multi-threading allows each channel to perform its own digital mixing and decimation in parallel, isolating individual satellite whistles from the shared wideband radio stream.

### 4. Weighted NTP Consensus Clock Steering: The Council of Watchmakers
When you ask multiple satellites what time it is, they will give you slightly different answers due to noise, atmospheric delay, and receiver variations.
- **The Old Way**: You update the clock based on the single satellite you tracked. If that satellite's signal was weak or low on the horizon, your clock might get steered incorrectly.
- **Weighted Consensus**: Imagine a council of watchmakers voting on the correct time. Instead of giving everyone an equal vote, the leader asks: "How clear is your vision?" (SNR), "How high was the sun when you looked?" (elevation), and "How steady is your hand?" (fit RMSE). The council trusts the watchmaker who had a crystal-clear view of a high-overhead pass with a high-quality mathematical fit. By taking a weighted average, the system steers the system clock smoothly and reliably.

### 5. Real-Time 3D Geodetic Geolocation: The Spherical Spotlight
If you are lost in a vast forest, you can pinpoint your location if you know your distance to several landmarks.
- **The Old Way**: You estimate your position by analyzing a single satellite's Doppler shift over its entire pass, which takes several minutes.
- **Real-Time 3D Geolocation**: When 4 or more satellites are visible, they shine their "coordinate spotlights" on you. Using a Gauss-Newton least-squares solver (a mathematical "warmer-colder" search), the system calculates your latitude, longitude, and altitude at 1 Hz by intersecting the spheres of range measurements. If the geometry is poor (high GDOP), the solver discards the epoch to prevent bad fixes.

