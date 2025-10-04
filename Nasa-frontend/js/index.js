
console.clear();

// Track variables
let scrollAnimation = null;
let currentVideo = "intro"; // start with intro.mp4
let videoElement = document.getElementById("backgroundVideo");
const container = document.querySelector("#container");

// Function to handle iOS play/pause
function handleIOSVideo() {
    videoElement.play().then(() => {
        videoElement.pause();
    }).catch(e => {
        console.log("iOS video handling:", e);
    });
}

// iOS specific handling
if (/iPad|iPhone|iPod/.test(navigator.userAgent)) {
    document.documentElement.addEventListener("touchstart", handleIOSVideo, { once: true });
}

// Initialize GSAP ScrollTrigger
gsap.registerPlugin(ScrollTrigger);

// Function to setup scroll animation for a video
function setupScrollAnimation(video) {
    if (scrollAnimation) scrollAnimation.kill();

    return gsap.timeline({
        scrollTrigger: {
            trigger: container,
            start: "top top",
            end: "bottom bottom",
            scrub: true,
            onUpdate: (self) => {
                const progress = self.progress;
                if (!isNaN(video.duration)) {
                    video.currentTime = video.duration * progress;
                }
            }
        }
    });
}

// Helper: switch video
function switchVideo(src, name) {
    if (currentVideo === name) return; // avoid duplicate switches

    // Save old video reference
    const oldVideo = videoElement;

    // Create new video element
    const newVideo = document.createElement("video");
    newVideo.id = "backgroundVideo";
    newVideo.className = "video-background";
    newVideo.src = src;
    newVideo.setAttribute("playsinline", "true");
    newVideo.muted = true;
    newVideo.preload = "auto";


    // Replace in DOM
    oldVideo.parentNode.replaceChild(newVideo, oldVideo);
    videoElement = newVideo; // update reference

    // Load new video
    newVideo.addEventListener('loadedmetadata', function () {
        newVideo.play().catch(e => console.log("Video play error:", e));
        scrollAnimation = setupScrollAnimation(newVideo);
    });

    // Handle error
    newVideo.addEventListener('error', function (e) {
        console.error("Video failed to load:", e);
        oldVideo.play();
    });

    currentVideo = name;
}

// Initial animation for intro video
videoElement.addEventListener('loadedmetadata', function () {
    scrollAnimation = setupScrollAnimation(videoElement);
}, { once: true });

// Unified scroll listener
window.addEventListener("scroll", () => {
    if (window.scrollY > 20) {
        switchVideo("../public/videos/gangster.mp4", "gangster");
    } else {
        switchVideo("../public/videos/intro.mp4", "intro");
    }

    // Navbar scroll effect
    const navbar = document.querySelector('.navbar');
    if (navbar) {
        if (window.scrollY > 50) {
            navbar.classList.add('scrolled');
        } else {
            navbar.classList.remove('scrolled');
        }
    }
});

// JavaScript for gaming elements
document.addEventListener('DOMContentLoaded', function () {
    // Nav menu interaction
    const navItems = document.querySelectorAll('.nav-item');
    navItems.forEach(item => {
        item.addEventListener('click', function (e) {
            e.preventDefault();
            navItems.forEach(i => i.classList.remove('active'));
            this.classList.add('active');
        });
    });

    // Toolbar interaction
    const toolbarItems = document.querySelectorAll('.toolbar-item');
    toolbarItems.forEach(item => {
        item.addEventListener('click', function () {
            toolbarItems.forEach(i => i.classList.remove('active'));
            this.classList.add('active');
        });
    });

    // Keyboard shortcuts for toolbar
    document.addEventListener('keydown', function (e) {
        const key = e.key.toUpperCase();
        const keyMap = {
            'R': 0, 'T': 1, 'I': 2, 'D': 3, 'A': 4, 'S': 5
        };

        if (key in keyMap) {
            e.preventDefault();
            const index = keyMap[key];
            if (toolbarItems[index]) {
                toolbarItems.forEach(i => i.classList.remove('active'));
                toolbarItems[index].classList.add('active');
            }
        }
    });
});
